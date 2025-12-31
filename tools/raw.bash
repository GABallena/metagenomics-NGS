#!/bin/bash
set -euo pipefail
shopt -s nullglob


# Directory containing raw reads
RAW_READS_DIR="raw_reads"  # Adjust this to your raw reads directory if different
FASTQC_OUTPUT_DIR="fastqc_output_raw"

# Activate the conda environment for FastQC
if command -v conda >/dev/null 2>&1; then
  # In some environments you may need: source "$HOME/miniconda3/etc/profile.d/conda.sh"
  conda activate fastqc_env >/dev/null 2>&1 || true
fi

# Create output directory if it doesn't exist
mkdir -p $FASTQC_OUTPUT_DIR

# Loop through each FASTQ file in the raw_reads directory
RAW_FILES=( "$RAW_READS_DIR"/*.fastq* )
if [ ${#RAW_FILES[@]} -eq 0 ]; then
  echo "No FASTQ files found in $RAW_READS_DIR"
  exit 0
fi

for RAW_FILE in $RAW_READS_DIR/SAMPLE*.fastq; do
  echo "Running FastQC on $RAW_FILE ..."
  
  # Run FastQC and direct the output to the specified directory
  fastqc -o $FASTQC_OUTPUT_DIR $RAW_FILE
  
  echo "FastQC complete for $RAW_FILE."
done

echo "All raw reads have been processed by FastQC."
