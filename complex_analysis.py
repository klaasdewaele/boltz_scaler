#!/usr/bin/env python3
"""
Complex Analysis Utility

A command line tool for analyzing protein complex predictions from Boltz.
Works with complexes of any size (monomers, dimers, trimers, tetramers, etc.).

Usage:
    python complex_analysis.py --fasta_files file1.fasta file2.fasta --predictions_dir /path/to/predictions 
                             [--complex_yamls_dir /path/to/yamls]
                             [--output_file results.csv]
                             [--plot_results]
"""

import os
import json
import yaml
import glob
import pandas as pd
import numpy as np
import re
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import logging
from Bio import SeqIO
import plotly.graph_objects as go
import plotly.io as pio
from plotly.offline import plot
import time
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("complex_analysis.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

class ComplexAnalyzer:
    def __init__(self, args):
        self.fasta_files = args.fasta_files
        self.predictions_dir = Path(args.predictions_dir)
        self.complex_yamls_dir = Path(args.complex_yamls_dir) if args.complex_yamls_dir else None
        self.output_file = args.output_file
        self.plot_results = args.plot_results

        # Check if paths exist
        if not all(os.path.exists(f) for f in self.fasta_files):
            missing = [f for f in self.fasta_files if not os.path.exists(f)]
            raise FileNotFoundError(f"FASTA file(s) not found: {', '.join(missing)}")
        
        if not self.predictions_dir.exists():
            raise FileNotFoundError(f"Predictions directory not found: {self.predictions_dir}")
        
        if self.complex_yamls_dir and not self.complex_yamls_dir.exists():
            raise FileNotFoundError(f"Complex YAML directory not found: {self.complex_yamls_dir}")
        
        # Storage for protein sequences
        self.sequence_to_gene = {}

    def extract_complex_number(self, path):
        """Extract the complex number from the path."""
        # Extract the base filename
        basename = os.path.basename(path)
        
        # First try to extract from the filename directly (most specific)
        # Looking for patterns like complex_123_model, tetramer_45_model, etc.
        match = re.search(r'(?:complex|tetramer|dimer|trimer|monomer|hexamer)_(\d+)_model', basename)
        if match:
            return match.group(1)
        
        # If that fails, extract from the directory name
        dirname = os.path.basename(os.path.dirname(path))
        match = re.search(r'(?:complex|tetramer|dimer|trimer|monomer|hexamer)_(\d+)', dirname)
        if match:
            return match.group(1)
        
        # Last resort: extract from the full path by finding the last occurrence of a number after any of the patterns
        matches = re.findall(r'(?:complex|tetramer|dimer|trimer|monomer|hexamer)_(\d+)', str(path))
        if matches:
            return matches[-1]  # Return the last match
        
        return None

    def load_sequences_from_fasta(self):
        """Load protein sequences from FASTA files and create a mapping from sequence to gene name."""
        logging.info("Loading sequences from FASTA files...")
        
        for fasta_file in self.fasta_files:
            logging.info(f"Processing FASTA file: {fasta_file}")
            try:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    gene_name = record.id
                    sequence = str(record.seq)
                    
                    # Clean up gene name (remove any pipe-separated parts, etc.)
                    if '|' in gene_name:
                        gene_name = gene_name.split('|')[0]
                    
                    gene_name = gene_name.strip()
                    self.sequence_to_gene[sequence] = gene_name
                    logging.info(f"Added sequence for gene: {gene_name}, length: {len(sequence)}")
            except Exception as e:
                logging.error(f"Error processing FASTA file {fasta_file}: {e}")
        
        logging.info(f"Total sequences loaded: {len(self.sequence_to_gene)}")
        return self.sequence_to_gene

    def find_yaml_file(self, complex_num):
        """Find the YAML file for a given complex number."""
        if not self.complex_yamls_dir:
            return None
        
        # Look for YAML files with various prefixes (complex_, tetramer_, etc.)
        patterns = [
            f"complex_{complex_num}.yaml", 
            f"tetramer_{complex_num}.yaml",
            f"dimer_{complex_num}.yaml",
            f"trimer_{complex_num}.yaml", 
            f"monomer_{complex_num}.yaml",
            f"hexamer_{complex_num}.yaml",
            # Add more patterns if needed
        ]
        
        for pattern in patterns:
            yaml_path = self.complex_yamls_dir / pattern
            if yaml_path.exists():
                return yaml_path
        
        # If specific files not found, try a glob pattern
        yaml_files = list(self.complex_yamls_dir.glob(f"*_{complex_num}.yaml"))
        if yaml_files:
            return yaml_files[0]
        
        return None

    def map_chains_to_genes(self, yaml_path):
        """Map chain IDs to gene names for a specific complex from a YAML file."""
        if not yaml_path or not os.path.exists(yaml_path):
            logging.warning(f"YAML file not found: {yaml_path}")
            return {}
        
        try:
            with open(yaml_path, 'r') as f:
                data = yaml.safe_load(f)
            
            chain_to_gene = {}
            
            if 'sequences' not in data:
                logging.warning(f"No sequences found in YAML file: {yaml_path}")
                return {}
            
            for protein in data['sequences']:
                if 'protein' not in protein:
                    continue
                
                chain_id = protein['protein']['id']
                sequence = protein['protein']['sequence']
                
                if sequence in self.sequence_to_gene:
                    gene_name = self.sequence_to_gene[sequence]
                    chain_to_gene[chain_id] = gene_name
                else:
                    # Try to extract gene name from the YAML if available
                    if 'id' in protein['protein'] and protein['protein']['id'] != chain_id:
                        gene_name = protein['protein']['id']
                        chain_to_gene[chain_id] = gene_name
                    else:
                        chain_to_gene[chain_id] = f"unknown_{chain_id}"
            
            return chain_to_gene
        
        except Exception as e:
            logging.error(f"Error processing YAML file {yaml_path}: {e}")
            return {}

    def extract_chain_ptms(self, json_data, chain_to_gene):
        """Extract pair-wise chain interface scores (pTMs) from JSON data."""
        try:
            chain_iptms = json_data.get('pair_chains_iptm', {})
            if not chain_iptms:
                return {}
            
            result_dict = {}
            
            # Create mapping from numeric indices to chain letters and gene names
            num_chains = len(chain_iptms)
            int_to_letter = {str(i): chr(65 + i) for i in range(num_chains)}  # 0->A, 1->B, etc.
            
            for chain_integer_str, interaction_dict in chain_iptms.items():
                chain_letter = int_to_letter[chain_integer_str]
                # Get the gene name for this chain
                if chain_letter not in chain_to_gene:
                    chain_name = f"unknown_{chain_letter}"
                else:
                    chain_name = chain_to_gene[chain_letter]
                
                for chain_integer_str_nested, iptm in interaction_dict.items():
                    if chain_integer_str_nested != chain_integer_str:  # Skip self-interactions
                        chain_letter_nested = int_to_letter[chain_integer_str_nested]
                        
                        # Get the gene name for the nested chain
                        if chain_letter_nested not in chain_to_gene:
                            chain_name_nested = f"unknown_{chain_letter_nested}"
                        else:
                            chain_name_nested = chain_to_gene[chain_letter_nested]
                        
                        # Create keys for this interaction
                        sorted_key = '_'.join(sorted([chain_name_nested, chain_name]))
                        unsorted_key = '_'.join([chain_name_nested, chain_name])
                        result_dict[unsorted_key] = iptm
            
            return result_dict
        
        except Exception as e:
            logging.error(f"Error extracting chain pTMs: {e}")
            return {}

    def process_complex_predictions(self):
        """Process all protein complex predictions and generate a DataFrame."""
        logging.info("Processing complex predictions...")
        
        # Prepare data for DataFrame
        data = []
        
        # Find all CIF files recursively
        pattern = f"{self.predictions_dir}/**/*_model_*.cif"
        cif_files = glob.glob(pattern, recursive=True)
        
        if not cif_files:
            logging.warning(f"No CIF files found with pattern: {pattern}")
            return pd.DataFrame()  # Return empty DataFrame
        
        logging.info(f"Found {len(cif_files)} CIF files to process")
        
        for cif_file in cif_files:
            # Extract the complex number
            complex_num = self.extract_complex_number(cif_file)
            if not complex_num:
                logging.warning(f"Could not extract complex number from path: {cif_file}")
                continue
            
            # Get the corresponding JSON confidence file
            cif_basename = os.path.basename(cif_file)
            json_path = os.path.join(os.path.dirname(cif_file), f"confidence_{cif_basename.replace('.cif', '.json')}")
            
            if not os.path.exists(json_path):
                logging.warning(f"No JSON confidence file found at: {json_path}")
                continue
            
            # Load the JSON data
            try:
                with open(json_path, 'r') as f:
                    json_data = json.load(f)
                
                # Find the YAML file for this complex
                yaml_path = self.find_yaml_file(complex_num)
                
                # Map chains to genes
                chain_to_gene = self.map_chains_to_genes(yaml_path)
                
                # Create a row for the DataFrame
                row_data = {
                    'complex_id': complex_num,
                    'confidence_score': json_data.get('confidence_score', None),
                    'ptm': json_data.get('ptm', None),
                    'iptm': json_data.get('iptm', None),
                    'complex_plddt': json_data.get('complex_plddt', None),
                }
                
                # Add chain-to-gene mapping
                for chain_id, gene_name in chain_to_gene.items():
                    row_data[f'chain_{chain_id}_gene'] = gene_name
                
                # Extract chain-pair pTMs
                chain_ptms = self.extract_chain_ptms(json_data, chain_to_gene)
                for key, value in chain_ptms.items():
                    row_data[key] = value
                
                # Add chain-specific PTM scores if available
                chains_ptm = json_data.get('chains_ptm', {})
                for chain_idx, ptm_value in chains_ptm.items():
                    # Try to convert numeric chain index to letter
                    try:
                        chain_letter = chr(65 + int(chain_idx))  # 0->A, 1->B, etc.
                        row_data[f'chain_{chain_letter}_ptm'] = ptm_value
                    except (ValueError, TypeError):
                        row_data[f'chain_{chain_idx}_ptm'] = ptm_value
                
                data.append(row_data)
            
            except Exception as e:
                logging.error(f"Error processing {json_path}: {e}")
        
        # Create DataFrame
        df = pd.DataFrame(data)
        
        # Sort by complex_id
        if not df.empty and 'complex_id' in df.columns:
            df['complex_id'] = df['complex_id'].astype(int)
            df = df.sort_values('complex_id')
        
        return df

    def create_plot(self, df):
        """Create a scatter plot visualizing the complex analysis results."""
        if df.empty:
            logging.warning("Cannot create plot: No data available")
            return
        
        logging.info("Creating visualization plot...")
        start_time = time.time()
        
        # Get all gene names from the dataframe
        gene_columns = [col for col in df.columns if col.startswith('chain_') and col.endswith('_gene')]
        all_genes = set()
        for col in gene_columns:
            all_genes.update(df[col].dropna().unique())
        
        # Create color mapping for genes (avoiding similar colors for adjacent genes)
        import colorsys
        
        # Generate distinct colors using HSV color space
        num_genes = len(all_genes)
        gene_colors = {}
        for i, gene in enumerate(sorted(all_genes)):
            if gene.startswith('unknown_'):  # Special case for unknown genes
                gene_colors[gene] = '#CCCCCC'  # Gray for unknown genes
            else:
                # Use HSV color space for better distribution
                hue = i / num_genes
                saturation = 0.7
                value = 0.85
                rgb = colorsys.hsv_to_rgb(hue, saturation, value)
                hex_color = "#{:02x}{:02x}{:02x}".format(
                    int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255)
                )
                gene_colors[gene] = hex_color
        
        # Find all gene pair columns
        gene_pair_columns = [col for col in df.columns if '_' in col and 
                           col not in gene_columns and 
                           not col.endswith('_ptm') and
                           not col.startswith('complex_')]
        
        # Data for plotting
        all_x = []  # Interface pTM (gene pair value)
        all_y = []  # Complex pTM
        all_colors = []  # Color for first gene
        all_border_colors = []  # Color for second gene
        all_gene_pairs = []  # Gene pair names
        all_customdata = []  # Data for hover info
        
        # Process data for plotting
        for _, row in df.iterrows():
            for col in gene_pair_columns:
                value = row[col]
                if pd.notna(value):
                    gene1, gene2 = col.split('_')
                    
                    # Add data
                    all_x.append(value)
                    all_y.append(row['ptm'])
                    all_colors.append(gene_colors.get(gene1, '#CCCCCC'))
                    all_border_colors.append(gene_colors.get(gene2, '#CCCCCC'))
                    all_gene_pairs.append(col)
                    
                    # Add custom data for hover info
                    custom_row = [
                        col,  # Gene pair
                        row['complex_id'],  # Complex ID
                        row.get('confidence_score', 'N/A'),  # Confidence score
                        row.get('iptm', 'N/A'),  # Interface pTM
                        row.get('complex_plddt', 'N/A')  # Complex pLDDT
                    ]
                    
                    # Add chain info to custom data
                    chains = []
                    for chain_col in gene_columns:
                        if pd.notna(row.get(chain_col)):
                            chains.append(str(row[chain_col]))
                    
                    custom_row.append(', '.join(chains))
                    all_customdata.append(custom_row)
        
        # Create figure
        fig = go.Figure()
        
        # Add all data points
        fig.add_trace(go.Scatter(
            x=all_x,
            y=all_y,
            mode='markers',
            marker=dict(
                color=all_colors,
                line=dict(
                    color=all_border_colors,
                    width=3
                ),
                size=8,
                opacity=0.7
            ),
            name='',
            showlegend=False,
            hovertemplate='<b>Gene Pair:</b> %{customdata[0]}<br>' +
                          '<b>Interface pTM:</b> %{x:.4f}<br>' +
                          '<b>Complex pTM:</b> %{y:.4f}<br>' +
                          '<b>Complex ID:</b> %{customdata[1]}<br>' +
                          '<b>Confidence Score:</b> %{customdata[2]:.4f}<br>' +
                          '<b>Interface pTM (iptm):</b> %{customdata[3]:.4f}<br>' +
                          '<b>Complex pLDDT:</b> %{customdata[4]:.4f}<br>' +
                          '<b>Chains:</b> %{customdata[5]}',
            customdata=all_customdata
        ))
        
        # Add legend entries for genes
        for gene, color in gene_colors.items():
            if gene.startswith('unknown_'):
                continue  # Skip unknown genes in legend
            
            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker=dict(
                    color=color,
                    size=10,
                    line=dict(width=0)
                ),
                name=gene
            ))
        
        # Update layout
        fig.update_layout(
            title='Complex Analysis Results',
            title_font=dict(size=20),
            xaxis=dict(
                title='Interface pTM',
                title_font=dict(size=16),
                range=[0, 1],
                gridcolor='lightgray',
                zeroline=True,
                zerolinecolor='lightgray'
            ),
            yaxis=dict(
                title='Complex pTM',
                title_font=dict(size=16),
                range=[0, 1],
                gridcolor='lightgray',
                zeroline=True,
                zerolinecolor='lightgray'
            ),
            plot_bgcolor='white',
            legend_title='Genes',
            legend_title_font=dict(size=14),
            hovermode='closest',
            width=1200,
            height=800,
            margin=dict(l=80, r=80, t=100, b=80)
        )
        
        # Add custom annotation
        fig.update_layout(
            annotations=[
                dict(
                    x=0.5,
                    y=-0.15,
                    xref="paper",
                    yref="paper",
                    text="Point fill color: First gene in pair | Point border color: Second gene in pair",
                    showarrow=False,
                    font=dict(size=12),
                    align="center"
                )
            ]
        )
        
        logging.info(f"Plot created in {time.time() - start_time:.2f} seconds")
        
        # Save the plot
        plot_filename = f"{os.path.splitext(self.output_file)[0]}_plot.html"
        plot(fig, filename=plot_filename, auto_open=False)
        logging.info(f"Plot saved to {plot_filename}")
        
        # Return the plot in case it's needed elsewhere
        return fig
    
    def run(self):
        """Execute the complete complex analysis pipeline."""
        # Step 1: Load sequences from FASTA files
        self.load_sequences_from_fasta()
        
        if not self.sequence_to_gene:
            logging.warning("No sequences loaded from FASTA files. Chain-to-gene mapping may be limited.")
        
        # Step 2: Process complex predictions
        df = self.process_complex_predictions()
        
        if df.empty:
            logging.error("No complex predictions were successfully processed.")
            return
        
        # Step 3: Sort by confidence score
        df_sorted = df.sort_values(by="confidence_score", ascending=False)
        
        # Step 4: Save results to CSV
        df_sorted.to_csv(self.output_file, index=False)
        logging.info(f"Results saved to {self.output_file}")
        
        # Step 5: Create visualization if requested
        if self.plot_results:
            self.create_plot(df_sorted)
        
        logging.info("Complex analysis completed successfully!")
        
        # Return the DataFrame for potential further analysis
        return df_sorted

def parse_args():
    parser = argparse.ArgumentParser(description="Complex Analysis Utility for Boltz Predictions")
    
    parser.add_argument("--fasta_files", nargs='+', required=True,
                        help="Path to one or more FASTA files containing protein sequences")
    
    parser.add_argument("--predictions_dir", required=True,
                        help="Directory containing the complex predictions")
    
    parser.add_argument("--complex_yamls_dir",
                        help="Directory containing the complex YAML files")
    
    parser.add_argument("--output_file", default="complex_analysis_results.csv",
                        help="Output CSV file for results (default: complex_analysis_results.csv)")
    
    parser.add_argument("--plot_results", action="store_true",
                        help="Generate an interactive visualization of the results")
    
    return parser.parse_args()

def main():
    """Main entry point for the complex analysis utility."""
    try:
        # Parse command line arguments
        args = parse_args()
        
        # Initialize and run the analyzer
        analyzer = ComplexAnalyzer(args)
        results = analyzer.run()
        
        # Print summary of results
        if results is not None and not results.empty:
            print("\nAnalysis Summary:")
            print(f"Total complexes analyzed: {len(results)}")
            print(f"Mean confidence score: {results['confidence_score'].mean():.4f}")
            print(f"Mean pTM score: {results['ptm'].mean():.4f}")
            if 'iptm' in results.columns:
                print(f"Mean iPTM score: {results['iptm'].mean():.4f}")
            
            # Print top 5 complexes by confidence score
            print("\nTop 5 complexes by confidence score:")
            top_5 = results.head(5)[['complex_id', 'confidence_score', 'ptm']]
            for _, row in top_5.iterrows():
                print(f"Complex {row['complex_id']}: Confidence={row['confidence_score']:.4f}, pTM={row['ptm']:.4f}")
            
            # Print gene breakdown
            chain_gene_columns = [col for col in results.columns if col.startswith('chain_') and col.endswith('_gene')]
            if chain_gene_columns:
                print("\nGenes found in complexes:")
                all_genes = set()
                for col in chain_gene_columns:
                    genes = results[col].dropna().unique()
                    all_genes.update(genes)
                
                for gene in sorted(all_genes):
                    if gene.startswith('unknown_'):
                        continue  # Skip unknown genes in summary
                    print(f"- {gene}")
        
        return 0
    
    except Exception as e:
        logging.error(f"An error occurred: {e}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())