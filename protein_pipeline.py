#!/Users/kdewaele/miniconda3/bin/python

"""
YAML-based Protein Structure Prediction Pipeline

This script takes FASTA files as input and creates a complete pipeline for:
1. Parsing the FASTA files and extracting unique protein sequences
2. Generating MSAs for each protein using YAML format
3. Creating YAML files for all protein combinations
4. Generating a batch script to run Boltz predictions

    This version uses YAML format exclusively to avoid FASTA parsing issues in Boltz.

    Usage:
        python yaml_protein_pipeline.py --fasta_files file1.fasta file2.fasta 
                                    [--r 2] 
                                    [--with_replacement] 
                                    [--skip_msa_generation]
                                    [--out_dir ./pipeline_output]
                                    [--recycling_steps 3]
                                    [--diffusion_samples 1]
                                    [--gpu_count 1]
                                    [--use_pbs]
                                    [--pbs_template path/to/template.pbs]
                                    [--pbs_work_dir path/to/work_dir]
                                    [--pbs_scripts_per_dir 1000]
                                    [--use_cpu]
                                    [--bait_proteins bait1.fasta bait2.fasta]
                                    [--candidate_proteins cand1.fasta cand2.fasta]
"""

import argparse
import itertools
import os
import re
import sys
import shutil
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional

try:
    from Bio import SeqIO
except ImportError:
    print("Biopython is required. Please install it with: pip install biopython")
    sys.exit(1)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("protein_pipeline.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

class YAMLProteinPipeline:
    def __init__(self, args):
        # Determine mode and set up file lists
        if args.mode == 'bait_candidate':
            self.mode = 'bait_candidate'
            self.bait_proteins = args.bait_proteins
            self.candidate_proteins = args.candidate_proteins
            self.fasta_files = list(set(args.bait_proteins + args.candidate_proteins))
        else:
            self.mode = 'all_vs_all'
            self.fasta_files = args.fasta_files
            self.bait_proteins = []
            self.candidate_proteins = []
        self.r = args.r
        self.with_replacement = args.with_replacement
        self.skip_msa_generation = args.skip_msa_generation
        self.out_dir = Path(args.out_dir)
        self.recycling_steps = args.recycling_steps
        self.diffusion_samples = args.diffusion_samples
        self.gpu_count = args.gpu_count
        self.use_cpu = args.use_cpu
        self.use_pbs = bool(args.use_pbs)
        if args.use_pbs and args.use_pbs is not True:
            self.pbs_template_path = args.use_pbs
        else:
            self.pbs_template_path = None
        self.pbs_work_dir = args.pbs_work_dir
        self.pbs_scripts_per_dir = args.pbs_scripts_per_dir
        
        # Set environment variable for CPU fallback if needed
        if self.use_cpu:
            os.environ['PYTORCH_ENABLE_MPS_FALLBACK'] = '1'
            logging.info("Setting PYTORCH_ENABLE_MPS_FALLBACK=1 for CPU usage")
        
        # Create necessary directories
        self.msas_dir = self.out_dir / "msas"
        self.yaml_dir = self.out_dir / "yaml_files"
        self.predictions_dir = self.out_dir / "predictions"
        self.pbs_dir = self.out_dir / "pbs_scripts"
        
        self.msas_dir.mkdir(parents=True, exist_ok=True)
        self.yaml_dir.mkdir(parents=True, exist_ok=True)
        self.predictions_dir.mkdir(parents=True, exist_ok=True)
        if self.use_pbs:
            self.pbs_dir.mkdir(parents=True, exist_ok=True)
        
        # Dict to store protein sequences and group info
        self.protein_info = {}
        
    def parse_fasta_files(self) -> dict:
        """Parse all FASTA files and extract unique protein sequences, tagging as bait/candidate if specified."""
        logging.info("Parsing FASTA files...")
        
        # Helper to tag group
        def get_group(fasta_file):
            if self.mode == 'bait_candidate':
                if fasta_file in self.bait_proteins:
                    return "bait"
                elif fasta_file in self.candidate_proteins:
                    return "candidate"
            return None
        
        for fasta_file in self.fasta_files:
            if not os.path.exists(fasta_file):
                logging.error(f"FASTA file not found: {fasta_file}")
                continue
            
            group = get_group(fasta_file)
            logging.info(f"Processing {fasta_file} (group: {group})")
            
            try:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    protein_id = record.id
                    sequence = str(record.seq)
                    
                    # If the ID contains a pipe character, use the part before it as the ID
                    if '|' in protein_id:
                        protein_id = protein_id.split('|')[0]
                    
                    # Ensure the ID is valid for file naming
                    protein_id = re.sub(r'[^\w]', '_', protein_id)
                    
                    # If this ID already exists with a different sequence, append a unique number
                    base_id = protein_id
                    counter = 1
                    while protein_id in self.protein_info and self.protein_info[protein_id]["sequence"] != sequence:
                        protein_id = f"{base_id}_{counter}"
                        counter += 1
                    
                    # Tag group if not already set
                    entry = self.protein_info.get(protein_id, {})
                    entry["sequence"] = sequence
                    if "group" not in entry or entry["group"] is None:
                        entry["group"] = group
                    self.protein_info[protein_id] = entry
                    logging.info(f"Added protein {protein_id} (group: {entry['group']}) with length {len(sequence)}")
            
            except Exception as e:
                logging.error(f"Error parsing {fasta_file}: {e}")
        
        logging.info(f"Total unique proteins found: {len(self.protein_info)}")
        return self.protein_info
    
    def generate_msas(self):
        """Generate MSAs for each protein using Boltz with YAML format."""
        if self.skip_msa_generation:
            logging.info("Skipping MSA generation as requested.")
            return
        
        logging.info("Generating MSAs for all proteins using YAML format...")
        
        for protein_id, sequence in self.protein_info.items():
            protein_dir = self.msas_dir / protein_id
            protein_dir.mkdir(exist_ok=True)
            
            # Use YAML format for better control
            yaml_path = protein_dir / f"{protein_id}.yaml"
            
            # Create a YAML file with simple, consistent chain IDs
            yaml_content = "version: 1\nsequences:\n"
            yaml_content += "  - protein:\n"
            yaml_content += "      id: A\n"  # Always use 'A' as the chain ID
            yaml_content += f"      sequence: {sequence}\n"
            
            with open(yaml_path, 'w') as f:
                f.write(yaml_content)
                
            logging.info(f"Generating MSA for {protein_id} using YAML format...")
            
            if self.use_cpu:
                logging.info("Using CPU for Boltz predictions")
                device_arg = "--accelerator cpu"
            else:
                logging.info("Using GPU for Boltz predictions")
                device_arg = f"--accelerator gpu --devices {self.gpu_count}"
            
            # Call Boltz to generate MSA
            cmd = [
                "boltz", "predict", 
                str(yaml_path), 
                "--use_msa_server", 
                "--out_dir", str(protein_dir / f"boltz_results_{protein_id}")
            ]
            
            # Add device argument
            cmd.extend(device_arg.split())
            
            try:
                subprocess.run(cmd, check=True)
                logging.info(f"MSA generation completed for {protein_id}")
            except subprocess.CalledProcessError as e:
                logging.error(f"MSA generation failed for {protein_id}: {e}")
                logging.error(f"Command used: {' '.join(cmd)}")
            except FileNotFoundError:
                logging.error("Boltz command not found. Please ensure Boltz is installed and in your PATH.")
                logging.info("Continuing with pipeline...")
    
    def scan_for_null_bytes(self, file_path):
        """Scan a file for null bytes and return their positions if found."""
        with open(file_path, 'rb') as f:
            content = f.read()
            null_positions = [i for i, byte in enumerate(content) if byte == 0]
            return null_positions, len(content)
    
    def fix_a3m_file(self, file_path, output_path=None):
        """
        Fix an A3M file by removing null bytes and ensuring it's properly formatted.
        
        Args:
            file_path: Path to the original A3M file
            output_path: Where to save the fixed file (defaults to original with '.fixed' suffix)
        
        Returns:
            dict: Summary of fixes applied
        """
        if output_path is None:
            output_path = str(file_path) + '.fixed'
        
        # Stats to track fixes
        fixes = {
            'null_bytes_removed': 0,
            'invalid_chars_removed': 0,
            'lines_processed': 0,
            'headers_processed': 0,
            'sequences_processed': 0
        }
        
        # Read as binary first to detect null bytes
        with open(file_path, 'rb') as f:
            binary_content = f.read()
            fixes['null_bytes_removed'] = binary_content.count(b'\x00')
            # Remove null bytes
            binary_content = binary_content.replace(b'\x00', b'')
        
        # Now process as text and fix other issues
        try:
            content = binary_content.decode('utf-8')
        except UnicodeDecodeError:
            # Fall back to latin-1 if UTF-8 fails
            content = binary_content.decode('latin-1')
        
        lines = content.split('\n')
        fixed_lines = []
        in_sequence = False
        
        for line in lines:
            fixes['lines_processed'] += 1
            stripped_line = line.strip()
            
            # Skip empty lines
            if not stripped_line:
                continue
            
            # Header line
            if stripped_line.startswith('>'):
                fixes['headers_processed'] += 1
                in_sequence = True
                fixed_lines.append(stripped_line)
                continue
            
            # Sequence line
            if in_sequence:
                fixes['sequences_processed'] += 1
                # Remove any characters that aren't valid amino acid codes or gaps
                # Allow both uppercase and lowercase letters (standard in A3M format)
                valid_sequence = re.sub(r'[^A-Za-z.-]', '', stripped_line)
                fixes['invalid_chars_removed'] += len(stripped_line) - len(valid_sequence)
                fixed_lines.append(valid_sequence)
        
        # Write the fixed content
        with open(output_path, 'w') as f:
            f.write('\n'.join(fixed_lines))
        
        return fixes
    
    def fix_msa_files(self):
        """Process all A3M files in the MSAs directory and fix any issues."""
        logging.info("Checking and fixing A3M files...")
        
        # Find all A3M files
        a3m_files = list(self.msas_dir.glob("**/*.a3m"))
        logging.info(f"Found {len(a3m_files)} A3M files to check")
        
        fixed_files = 0
        clean_files = 0
        total_null_bytes = 0
        
        # Process each file
        for file_path in a3m_files:
            logging.info(f"Checking {file_path}")
            null_positions, file_size = self.scan_for_null_bytes(file_path)
            
            if null_positions:
                fixed_path = str(file_path) + '.fixed'
                logging.warning(f"Found {len(null_positions)} null bytes in {file_path}")
                total_null_bytes += len(null_positions)
                
                # Create a backup
                backup_path = str(file_path) + '.bak'
                shutil.copy2(file_path, backup_path)
                logging.info(f"Created backup at {backup_path}")
                
                # Fix the file
                fixes = self.fix_a3m_file(file_path, fixed_path)
                logging.info(f"Fixed file saved to {fixed_path}")
                
                # Replace the original with the fixed version
                shutil.move(fixed_path, file_path)
                logging.info(f"Replaced original with fixed version")
                
                fixed_files += 1
                logging.info(f"Fixes applied: {fixes}")
            else:
                clean_files += 1
                logging.info(f"No null bytes found, file appears clean")
        
        # Print summary
        logging.info("========== MSA FIX SUMMARY ==========")
        logging.info(f"Total A3M files checked: {len(a3m_files)}")
        logging.info(f"Clean files: {clean_files}")
        logging.info(f"Fixed files: {fixed_files}")
        logging.info(f"Total null bytes removed: {total_null_bytes}")
    
    def get_msa_path(self, protein_id):
        """Get the appropriate MSA file path for a protein, adjusted for YAML approach."""
        # When using YAML input with chain ID 'A', the MSA paths are different
        msa_path = self.msas_dir / protein_id / f"boltz_results_{protein_id}" / f"boltz_results_{protein_id}" / "msa" / f"{protein_id}_unpaired_tmp_env" / "bfd.mgnify30.metaeuk30.smag30.a3m"

        # Check if file exists
        if msa_path.exists():
            path_str = str(msa_path)
            if self.use_pbs:
                # Map local path under out_dir to HPC path under pbs_work_dir
                rel = os.path.relpath(path_str, start=str(self.out_dir))
                return os.path.join(self.pbs_work_dir, rel)
            return path_str
        else:
            # Fallback to uniref if BFD doesn't exist
            logging.info(f'Falling back to uniref for {protein_id}')
            backup_path = self.msas_dir / protein_id / f"boltz_results_{protein_id}" / f"boltz_results_{protein_id}" / "msa" / "A_unpaired_tmp_env" / "uniref.a3m"
            if backup_path.exists():
                logging.warning(f"Using backup MSA for {protein_id}")
                path_str = str(backup_path)
                if self.use_pbs:
                    rel = os.path.relpath(path_str, start=str(self.out_dir))
                    return os.path.join(self.pbs_work_dir, rel)
                return path_str
            else:
                logging.warning(f"No MSA found for {protein_id}")
                return "empty"  # Use 'empty' as fallback if no MSA found
    
    def generate_yaml_files(self):
        """Generate YAML files for protein combinations, supporting bait/candidate logic if specified."""
        logging.info("Generating YAML files for protein combinations...")
        
        # Get all protein IDs and their groups
        protein_ids = list(self.protein_info.keys())
        bait_ids = [pid for pid in protein_ids if self.protein_info[pid].get("group") == "bait"]
        candidate_ids = [pid for pid in protein_ids if self.protein_info[pid].get("group") == "candidate"]
        
        # If bait/candidate mode is active and both lists are non-empty, use special logic
        if bait_ids and candidate_ids:
            logging.info(f"Bait/candidate mode: {len(bait_ids)} baits, {len(candidate_ids)} candidates")
            combinations = []
            # Bait-candidate combinations
            for bait in bait_ids:
                for candidate in candidate_ids:
                    if self.r == 2:
                        combinations.append((bait, candidate))
                    else:
                        # For r > 2, fill with baits and one candidate
                        for bait_combo in itertools.combinations([b for b in bait_ids if b != bait], self.r - 2):
                            combinations.append((bait, candidate) + bait_combo)
        else:
            # Default: all-vs-all logic
            if self.with_replacement:
                combinations = list(itertools.combinations_with_replacement(protein_ids, self.r))
            else:
                combinations = list(itertools.combinations(protein_ids, self.r))
        
        logging.info(f"Total combinations: {len(combinations)}")
        
        # Generate YAML for each combination
        for i, combo in enumerate(combinations):
            yaml_content = "version: 1\nsequences:\n"
            
            # Add each protein in the combination
            for j, protein_id in enumerate(combo):
                chain_id = chr(65 + j)  # A, B, C, D, etc.
                yaml_content += f"  - protein:\n"
                yaml_content += f"      id: {chain_id}\n"
                yaml_content += f"      sequence: {self.protein_info[protein_id]['sequence']}\n"
                yaml_content += f"      msa: {self.get_msa_path(protein_id)}\n"
            
            # Write the YAML file
            yaml_file = self.yaml_dir / f"complex_{i+1}.yaml"
            with open(yaml_file, "w") as f:
                f.write(yaml_content)
        
        logging.info(f"Generated {len(combinations)} YAML files in {self.yaml_dir}")
    
    def create_batch_script(self):
        """Create a batch processing script for Boltz predictions."""
        logging.info("Creating batch processing script...")
        
        batch_script = f"""#!/bin/bash
# Batch processing script for Boltz predictions

# Directory containing all YAML files
YAML_DIR="{self.yaml_dir}"
OUT_DIR="{self.predictions_dir}"

# Create output directory
mkdir -p $OUT_DIR

# Get total number of files
TOTAL=$(ls $YAML_DIR/*.yaml | wc -l)
echo "Found $TOTAL YAML files to process"

# Process options
RECYCLING="{self.recycling_steps}"  # Recycling steps
SAMPLES="{self.diffusion_samples}"  # Diffusion samples
GPU_COUNT="{self.gpu_count}"        # Number of GPUs to use

echo "Using parameters: --recycling_steps $RECYCLING --diffusion_samples $SAMPLES --devices $GPU_COUNT"

# Process all files
for yaml_file in $YAML_DIR/*.yaml; do
    base_name=$(basename "$yaml_file" .yaml)
    echo "Processing $base_name..."
    
    # Run Boltz prediction
    if [ -d "$OUT_DIR/$base_name" ] && [ ! -z "$(ls -A $OUT_DIR/$base_name)" ]; then
        echo "Output already exists for $base_name, skipping"
    else
        boltz predict "$yaml_file" --out_dir "$OUT_DIR/$base_name" --recycling_steps $RECYCLING --diffusion_samples $SAMPLES --devices $GPU_COUNT
    fi
done

echo "All predictions completed!"
"""
        
        # Write the batch script
        script_path = self.out_dir / "run_predictions.sh"
        with open(script_path, "w") as f:
            f.write(batch_script)
        
        # Make the script executable
        os.chmod(script_path, 0o755)
        
        logging.info(f"Created batch processing script '{script_path}'")
    
    def read_pbs_template(self) -> str:
        """Read and validate the PBS/SLURM template file."""
        import glob
        if self.pbs_template_path:
            template_path = Path(self.pbs_template_path)
            if not template_path.exists():
                raise FileNotFoundError(f"PBS/SLURM template file not found: {template_path}")
        else:
            # Search out_dir for a file containing 'pbs_template' or 'slurm_template'
            candidates = list(self.out_dir.glob("*pbs_template*")) + list(self.out_dir.glob("*slurm_template*"))
            if not candidates:
                raise FileNotFoundError("No PBS/SLURM template file found in output directory. Please provide a template file or specify it as an argument to --use_pbs.")
            template_path = candidates[0]
        with open(template_path) as f:
            template = f.read()
        # Validate that template has required PBS/SLURM directives
        if '#PBS' not in template and '#SBATCH' not in template:
            raise ValueError("Template file does not appear to be a valid PBS or SLURM script.")
        return template
    
    def create_pbs_scripts(self):
        """Create individual PBS scripts for each complex prediction."""
        logging.info("Creating PBS scripts for complex predictions...")
        
        # Read the template
        template = self.read_pbs_template()
        
        # Get list of all YAML files
        yaml_files = list(self.yaml_dir.glob("*.yaml"))
        total_files = len(yaml_files)
        
        # Calculate number of subdirectories needed
        num_subdirs = (total_files + self.pbs_scripts_per_dir - 1) // self.pbs_scripts_per_dir
        
        logging.info(f"Creating PBS scripts for {total_files} complex predictions")
        logging.info(f"Splitting into {num_subdirs} subdirectories")
        
        for i, yaml_file in enumerate(yaml_files):
            # Determine subdirectory
            subdir_num = i // self.pbs_scripts_per_dir
            pbs_subdir = self.pbs_dir / f"batch_{subdir_num + 1}"
            pbs_subdir.mkdir(exist_ok=True)
            
            complex_name = yaml_file.stem
            complex_num = int(complex_name.split('_')[1])  # Extract number from "complex_1"
            
            # Create prediction output directory path for HPC environment
            hpc_out_dir = f"{self.pbs_work_dir}/predictions/{complex_name}"

            # Create cache directory path for HPC environment
            hpc_cache_dir = f"{self.pbs_work_dir}/cache"
            
            # Replace placeholders in template
            script_content = template.replace('{complex_num}', str(complex_num))
            
            # Add our custom commands after the last PBS directive
            device_arg = "--accelerator cpu" if self.use_cpu else f"--devices {self.gpu_count}"
            
            custom_commands = f"""
# Create output directory
mkdir -p {hpc_out_dir}

# Set environment variable for CPU fallback if needed
{f'export PYTORCH_ENABLE_MPS_FALLBACK=1' if self.use_cpu else ''}

# Make sure MSA and YAML files are available on HPC server

# Run Boltz prediction
boltz predict "{self.pbs_work_dir}/yaml_files/{yaml_file.name}" \\
    --out_dir "{hpc_out_dir}" \\
    --recycling_steps {self.recycling_steps} \\
    --diffusion_samples {self.diffusion_samples} \\
    {device_arg} \\
    --cache {hpc_cache_dir}
"""
            
            # Find the last line
            lines = script_content.split('\n')
            
            # Insert our commands after the last line of template file
            lines.insert(len(lines), custom_commands)
            script_content = '\n'.join(lines)
            
            # Write the PBS script
            script_path = pbs_subdir / f"{complex_name}.pbs"
            with open(script_path, 'w') as f:
                f.write(script_content)
            
            # Make script executable
            os.chmod(script_path, 0o755)
        
        logging.info(f"Created {total_files} PBS scripts in {self.pbs_dir}")
        logging.info("To submit jobs, use: qsub path/to/script.pbs")
    
    def run(self):
        """Execute the complete pipeline."""
        logging.info("Starting YAML-based protein structure prediction pipeline...")
        
        # Step 1: Parse FASTA files
        self.parse_fasta_files()
        
        # Step 2: Generate MSAs
        self.generate_msas()
        
        # Step 3: Fix any issues in MSA files
        self.fix_msa_files()
        
        # Step 4: Generate YAML files
        self.generate_yaml_files()
        
        # Step 5: Create either batch script or PBS scripts
        if self.use_pbs:
            self.create_pbs_scripts()
            logging.info("Pipeline completed successfully!")
            logging.info(f"PBS scripts directory: {self.pbs_dir}")
            logging.info("Before submitting jobs:")
            logging.info(f"1. Copy integral working directory to HPC environment (comprising dirs: msas, pbs_scripts, predictions, yaml_files")
            logging.info("2. Submit jobs using: for f in path/to/pbs_scripts/*/*.pbs; do qsub $f; done")
        else:
            self.create_batch_script()
            logging.info("Pipeline completed successfully!")
            logging.info(f"Output directory: {self.out_dir}")
            logging.info(f"To run predictions, execute: {self.out_dir}/run_predictions.sh")


def parse_args():
    parser = argparse.ArgumentParser(description="YAML-based Protein Structure Prediction Pipeline")
    
    parser.add_argument("--fasta_files", nargs='+',
                        help="Path to one or more FASTA files for all-vs-all combinations")
    parser.add_argument("--bait_proteins", nargs='+',
                        help="List of bait protein FASTA files (for bait/candidate mode)")
    parser.add_argument("--candidate_proteins", nargs='+',
                        help="List of candidate protein FASTA files (for bait/candidate mode)")
    parser.add_argument("--r", type=int, default=2,
                        help="Number of elements in each combination (default: 2)")
    parser.add_argument("--with_replacement", action="store_true",
                        help="Use combinations with replacement (default: without replacement)")
    parser.add_argument("--skip_msa_generation", action="store_true",
                        help="Skip MSA generation (useful if MSAs already exist)")
    parser.add_argument("--out_dir", default="./yaml_pipeline_output",
                        help="Output directory for all generated files (default: ./yaml_pipeline_output)")
    parser.add_argument("--recycling_steps", type=int, default=3,
                        help="Number of recycling steps for Boltz predictions (default: 3)")
    parser.add_argument("--diffusion_samples", type=int, default=1,
                        help="Number of diffusion samples for Boltz predictions (default: 1)")
    parser.add_argument("--gpu_count", type=int, default=1,
                        help="Number of GPUs to use for Boltz predictions (default: 1)")
    parser.add_argument("--use_pbs", nargs="?", const=True, default=False,
                        help="Generate PBS scripts instead of a single batch script. Optionally provide a template path.")
    parser.add_argument("--pbs_work_dir", type=str, default="~/boltz_predictions",
                        help="Working directory on HPC system (default: home directory)")
    parser.add_argument("--pbs_scripts_per_dir", type=int, default=1000,
                        help="Maximum number of PBS scripts per directory (default: 1000)")
    parser.add_argument("--use_cpu", action="store_true",
                        help="Use CPU instead of GPU for Boltz predictions")
    
    args = parser.parse_args()
    
    # Infer mode and validate
    if args.bait_proteins and args.candidate_proteins:
        if args.fasta_files:
            print("Warning: --fasta_files is ignored when both --bait_proteins and --candidate_proteins are provided.")
        mode = 'bait_candidate'
    elif args.fasta_files:
        mode = 'all_vs_all'
    else:
        parser.error("You must provide either --fasta_files for all-vs-all mode, or both --bait_proteins and --candidate_proteins for bait/candidate mode.")
    args.mode = mode
    return args


if __name__ == "__main__":
    args = parse_args()
    pipeline = YAMLProteinPipeline(args)
    pipeline.run()