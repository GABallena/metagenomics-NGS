#!/bin/bash

# Strict mode
set -euo pipefail
IFS=$'\n\t'

# Function to display usage
usage() {
    cat << EOF
Usage: $(basename "$0") MERGED_READS READ1 READ2 CONTIGS OUTPUT_DIR

Process and filter sequencing reads and contigs.

Arguments:
    MERGED_READS  Path to merged reads (FASTQ)
    READ1         Path to original R1 reads
    READ2         Path to original R2 reads
    CONTIGS       Path to contigs file (FASTA)
    OUTPUT_DIR    Directory for outputs
EOF
    exit 1
}

# Function to check if file exists and is not empty
check_file() {
    if [[ ! -f "$1" || ! -s "$1" ]]; then
        echo "Error: File '$1' does not exist or is empty" >&2
        exit 1
    fi
}

# Function to process files
process_files() {
    local merged_reads=$1
    local read1=$2
    local read2=$3
    local contigs=$4
    local output_dir=$5
    local log_dir="$output_dir/logs"

    # Define output files
    local paired_r1="$output_dir/paired_reads_paired_R1.fastq"
    local paired_r2="$output_dir/paired_reads_paired_R2.fastq"
    local processed_contigs="$output_dir/processed_contigs.fa"
    local split_r1="$output_dir/longreads_R1.fastq"
    local split_r2="$output_dir/longreads_R2.fastq"
    local merged_headers="$log_dir/merged_reads.headers"

    # Create output directories
    mkdir -p "$output_dir" "$log_dir"

    echo "Extracting headers from merged reads..."
    grep "^@" "$merged_reads" | cut -f1 -d' ' | sed 's/^@//' > "$merged_headers"

    echo "Filtering R1 and R2 based on headers..."
    seqtk subseq "$read1" "$merged_headers" > "$paired_r1"
    seqtk subseq "$read2" "$merged_headers" > "$paired_r2"
    check_file "$paired_r1"
    check_file "$paired_r2"

    echo "Reformatting contigs..."
    seqtk seq -l 0 "$contigs" > "$processed_contigs"
    check_file "$processed_contigs"

    echo "Splitting merged reads into R1 and R2..."
    seqtk seq -l0 -1 "$merged_reads" > "$split_r1"
    seqtk seq -l0 -2 "$merged_reads" > "$split_r2"
    check_file "$split_r1"
    check_file "$split_r2"
}

# Check number of arguments
[[ $# -eq 5 ]] || usage

# Check if input files exist
for file in "$1" "$2" "$3" "$4"; do
    check_file "$file"
done

# Check if seqtk is available
if ! command -v seqtk >/dev/null 2>&1; then
    echo "Error: seqtk is required but not found in PATH" >&2
    exit 1
fi

# Process files
process_files "$1" "$2" "$3" "$4" "$5"

echo "All tasks completed successfully."
