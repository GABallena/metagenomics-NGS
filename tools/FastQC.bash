#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Usage: bash FASTQC.bash <input_R1> <input_R2> <output_dir>

# Ensure the script exits on errors or unset variables
set -e
set -u

# Input files and output directory
INPUT_R1=$1
INPUT_R2=$2
OUTPUT_DIR=$3

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Extract the base name (e.g., SSR1, SSR2) from the input files
SAMPLE_ID="$(basename "$INPUT_R1")"
SAMPLE_ID="${SAMPLE_ID%%_R1.fastq.gz}"
SAMPLE_ID="${SAMPLE_ID%%_1.fastq.gz}"

echo "Running FastQC for $SAMPLE_ID ..."

# Run FastQC on both paired-end files
fastqc -o "$OUTPUT_DIR" "$INPUT_R1" "$INPUT_R2"

echo "FastQC for $SAMPLE_ID complete."

echo "All FastQC outputs for $SAMPLE_ID are stored in $OUTPUT_DIR."
