#!/usr/bin/env python3
"""
Script to align multiple FASTA files and generate a phylogenetic tree using MAFFT, TrimAl, and FastTree.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

def run_command(cmd, description):
    """Run a shell command and handle errors."""
    print(f"Running: {description}")
    print(f"Command: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=True, 
                              capture_output=True, text=True)
        print("Success!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running {description}:")
        print(f"Return code: {e.returncode}")
        print(f"Stderr: {e.stderr}")
        print(f"Stdout: {e.stdout}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Align FASTA files and generate phylogenetic trees using MAFFT, TrimAl, and FastTree"
    )
    
    parser.add_argument(
        "-i", "--input", 
        nargs="+", 
        required=True,
        help="Input FASTA files (one or more)"
    )
    
    parser.add_argument(
        "-o", "--outdir", 
        required=True,
        help="Output directory"
    )
    
    parser.add_argument(
        "-n", "--name", 
        default="merged_selfasttree_Flavodoxin",
        help="Base name for output files (default: merged_selfasttree_Flavodoxin)"
    )
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    
    # Check if input files exist
    for fasta_file in args.input:
        if not os.path.exists(fasta_file):
            print(f"Error: Input file {fasta_file} does not exist!")
            sys.exit(1)
    
    # Define file paths
    FASTA = os.path.join(args.outdir, f"{args.name}.fasta")
    ALIGNED = FASTA.replace(".fasta", "_align.fasta")
    TRIM = ALIGNED.replace("_align.fasta", "_trimed.fasta")
    FASTTREE = ALIGNED.replace("_align.fasta", "_fasttree.nw")
    
    print(f"Input files: {args.input}")
    print(f"Output directory: {args.outdir}")
    print(f"Base name: {args.name}")
    print(f"FASTA: {FASTA}")
    print(f"ALIGNED: {ALIGNED}")
    print(f"TRIM: {TRIM}")
    print(f"FASTTREE: {FASTTREE}")
    
    # Step 1: Merge input FASTA files (if multiple files provided)
    if len(args.input) > 1:
        print("\nMerging FASTA files...")
        with open(FASTA, 'w') as outfile:
            for fasta_file in args.input:
                with open(fasta_file, 'r') as infile:
                    outfile.write(infile.read())
    else:
        # If only one file, create a symlink or copy
        os.symlink(os.path.abspath(args.input[0]), FASTA)
    
    # Step 2: Run MAFFT alignment
    mafft_cmd = f"mafft {FASTA} > {ALIGNED}"
    if not run_command(mafft_cmd, "MAFFT alignment"):
        sys.exit(1)
    
    # Step 3: Run TrimAl
    trimal_cmd = f"trimal -in {ALIGNED} -out {TRIM} -gt 0.3"
    if not run_command(trimal_cmd, "TrimAl trimming"):
        sys.exit(1)
    
    # Step 4: Run FastTree
    fasttree_cmd = f"fasttree {TRIM} > {FASTTREE}"
    if not run_command(fasttree_cmd, "FastTree phylogeny"):
        sys.exit(1)
    
    print(f"\nPipeline completed successfully!")
    print(f"Alignment: {ALIGNED}")
    print(f"Trimmed alignment: {TRIM}")
    print(f"Phylogenetic tree: {FASTTREE}")

if __name__ == "__main__":
    main()
