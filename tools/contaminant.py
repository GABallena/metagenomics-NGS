#!/usr/bin/env python3
import argparse
import subprocess
import os
import shutil

def parse_args():
    parser = argparse.ArgumentParser(description='Remove contaminant reads using bowtie2')
    parser.add_argument('-1', '--input1', required=True, help='Input fastq file 1')
    parser.add_argument('-2', '--input2', required=True, help='Input fastq file 2')
    parser.add_argument('-x', '--index', required=True, help='Bowtie2 index prefix')
    parser.add_argument('-o1', '--output1', required=True, help='Output clean fastq file 1')
    parser.add_argument('-o2', '--output2', required=True, help='Output clean fastq file 2')
    parser.add_argument('-s', '--sam', required=True, help='Output SAM file')
    parser.add_argument('-t', '--threads', default=1, type=int, help='Number of threads')
    parser.add_argument('--input_files', nargs='+', help='List of input FASTA files to concatenate')
    parser.add_argument('--merged_output', help='Path for merged FASTA file')
    parser.add_argument('--db_prefix', help='Prefix for bowtie2 database')
    return parser.parse_args()

def prepare_contaminant_db(input_files, merged_output, db_prefix):
    """Prepare contaminant database by merging files and building bowtie2 index"""
    # Create output directory
    os.makedirs(os.path.dirname(merged_output), exist_ok=True)
    
    # Concatenate files
    with open(merged_output, 'w') as outfile:
        for fname in input_files:
            with open(fname) as infile:
                shutil.copyfileobj(infile, outfile)
    
    # Build bowtie2 index
    cmd = ['bowtie2-build', '--quiet', merged_output, db_prefix]
    subprocess.run(cmd, check=True)

def map_and_filter(args):
    """Map reads and filter contaminants"""
    # Create output directory
    os.makedirs(os.path.dirname(args.output1), exist_ok=True)
    
    # Map reads to contaminant database using Bowtie2
    bowtie2_cmd = [
        "bowtie2",
        "-p", str(args.threads),
        "-x", args.index,
        "-1", args.input1,
        "-2", args.input2,
        "-S", args.sam,
        "--un-conc", os.path.join(os.path.dirname(args.output1), "temp_clean_" + str(os.getpid()))
    ]
    
    subprocess.run(bowtie2_cmd, check=True)
    
    # Rename the output files
    os.rename(os.path.join(os.path.dirname(args.output1), "temp_clean.1"), args.output1)
    os.rename(os.path.join(os.path.dirname(args.output1), "temp_clean.2"), args.output2)

def main():
    args = parse_args()
    
    # If preparing contaminant database
    if args.input_files and args.merged_output and args.db_prefix:
        try:
            prepare_contaminant_db(args.input_files, args.merged_output, args.db_prefix)
            print(f"Successfully created contaminant database at {args.db_prefix}")
            return
        except Exception as e:
            print(f"Error creating contaminant database: {str(e)}")
            raise
    
    # Check if input files exist
    for f in [args.input1, args.input2]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input file {f} not found")
    
    # Run contaminant removal
    map_and_filter(args)
    
    print("Contaminant removal complete")
    print(f"Clean reads written to {args.output1} and {args.output2}")
    print(f"Mapping results written to {args.sam}")

if __name__ == "__main__":
    main()