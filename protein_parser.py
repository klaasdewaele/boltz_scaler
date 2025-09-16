#!/usr/bin/env python3
"""
Shared protein parsing utilities for Boltz pipeline

This module provides consistent FASTA/text file parsing logic that is used by both
the command-line tools and the GUI to ensure maintainability and consistency.
"""

import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, NamedTuple

try:
    from Bio import SeqIO
except ImportError:
    print("Biopython is required. Please install it with: pip install biopython")
    raise


class ProteinRecord(NamedTuple):
    """Data structure to represent a protein record"""
    id: str
    name: str  # Original record ID from file
    sequence: str
    source_file: str
    group: Optional[str] = None  # For bait/candidate grouping


def parse_protein_files(file_paths: List[str], mode: str = 'standard') -> Dict[str, ProteinRecord]:
    """
    Parse protein files (FASTA or text) and extract unique protein sequences.
    
    This is the shared parsing logic used by both CLI tools and GUI.
    
    Args:
        file_paths: List of file paths to parse
        mode: Parsing mode ('standard' or 'bait_candidate')
        
    Returns:
        Dictionary mapping unique protein IDs to ProteinRecord objects
    """
    protein_info = {}
    
    for file_path in file_paths:
        if not os.path.exists(file_path):
            logging.warning(f"File not found: {file_path}")
            continue
            
        logging.info(f"Processing {file_path}")
        
        try:
            # Determine file format
            file_ext = Path(file_path).suffix.lower()
            
            if file_ext == '.txt':
                records = _parse_text_file(file_path)
            else:
                # Handle FASTA files (.fasta, .fas, .fa)
                records = _parse_fasta_file(file_path)
            
            # Process records and ensure unique IDs
            for original_id, sequence in records:
                # Clean up protein ID (same logic as CLI tools)
                protein_id = original_id
                
                # If the ID contains a pipe character, use the part before it as the ID
                if '|' in protein_id:
                    protein_id = protein_id.split('|')[0]
                
                # Ensure the ID is valid for file naming
                protein_id = re.sub(r'[^\w]', '_', protein_id)
                
                # If this ID already exists with a different sequence, append a unique number
                base_id = protein_id
                counter = 1
                while protein_id in protein_info and protein_info[protein_id].sequence != sequence:
                    protein_id = f"{base_id}_{counter}"
                    counter += 1
                
                # Create protein record
                protein_record = ProteinRecord(
                    id=protein_id,
                    name=original_id,
                    sequence=sequence,
                    source_file=file_path
                )
                
                protein_info[protein_id] = protein_record
                logging.debug(f"Added protein {protein_id} from {file_path}")
                
        except Exception as e:
            logging.error(f"Error processing {file_path}: {str(e)}")
            continue
    
    logging.info(f"Successfully parsed {len(protein_info)} unique proteins from {len(file_paths)} files")
    return protein_info


def _parse_fasta_file(file_path: str) -> List[Tuple[str, str]]:
    """Parse a FASTA file and return list of (id, sequence) tuples"""
    records = []
    
    with open(file_path, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            records.append((record.id, str(record.seq)))
    
    return records


def _parse_text_file(file_path: str) -> List[Tuple[str, str]]:
    """Parse a text file and return list of (id, sequence) tuples"""
    records = []
    
    with open(file_path, 'r') as f:
        content = f.read().strip()
        
        if not content:
            return records
        
        # Try to parse as FASTA first
        try:
            f.seek(0)
            fasta_records = list(SeqIO.parse(f, "fasta"))
            if fasta_records:
                # It's FASTA format in a txt file
                for record in fasta_records:
                    records.append((record.id, str(record.seq)))
            else:
                # Plain text format - treat as single sequence
                file_stem = Path(file_path).stem
                clean_sequence = content.replace('\n', '').replace(' ', '').replace('\t', '')
                records.append((f"protein_from_{file_stem}", clean_sequence))
                
        except Exception:
            # Plain text format - treat as single sequence
            file_stem = Path(file_path).stem
            clean_sequence = content.replace('\n', '').replace(' ', '').replace('\t', '')
            records.append((f"protein_from_{file_stem}", clean_sequence))
    
    return records


def get_protein_combinations(proteins: Dict[str, ProteinRecord], 
                           r: int, 
                           with_replacement: bool = False) -> List[Tuple[str, ...]]:
    """
    Generate protein combinations for complex prediction.
    
    Args:
        proteins: Dictionary of protein records
        r: Number of proteins per combination
        with_replacement: Allow repeated proteins in combinations
        
    Returns:
        List of protein ID tuples representing combinations
    """
    import itertools
    
    protein_ids = list(proteins.keys())
    
    if with_replacement:
        combinations = list(itertools.combinations_with_replacement(protein_ids, r))
    else:
        combinations = list(itertools.combinations(protein_ids, r))
    
    return combinations


# For backward compatibility with existing code
def parse_fasta_files_legacy(file_paths: List[str]) -> Dict[str, str]:
    """
    Legacy function that returns the old format (protein_id -> sequence mapping).
    Used for backward compatibility with existing CLI tools.
    """
    protein_records = parse_protein_files(file_paths)
    return {pid: record.sequence for pid, record in protein_records.items()}


if __name__ == "__main__":
    # Simple test
    import sys
    if len(sys.argv) > 1:
        test_files = sys.argv[1:]
        logging.basicConfig(level=logging.INFO)
        
        proteins = parse_protein_files(test_files)
        print(f"\nFound {len(proteins)} proteins:")
        for pid, record in proteins.items():
            print(f"  {pid}: {record.name} ({len(record.sequence)} aa)")