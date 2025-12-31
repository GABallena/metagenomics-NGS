#!/bin/bash
set -euo pipefail

# Public-safe template:
# Set these to your local database FASTA paths (not provided in this repository).
NCBI_AMR_FASTA="${NCBI_AMR_FASTA:-databases/amr_db_a.fasta}"
CARD_FASTA="${CARD_FASTA:-databases/amr_db_b.fasta}"


# Define the directories
CONTIGS_DIR="merged_contigs"      # Directory with all the filtered contigs
OUTPUT_DIR="prokka_output/prokka_ARGs"   # Directory for Prokka output
EXTRA_DIR="databases/arg_proteins.fasta" # Protein FASTA file used by Prokka (template path)

# Concatenate all ARG nucleotide sequences into a single file
cat "$NCBI_AMR_FASTA" \
    "$CARD_FASTA" > all_ARG_sequences.fasta

# Translate the concatenated nucleotide sequences into proteins
# If your input FASTA is nucleotide sequences, translate to proteins:
# transeq -sequence all_ARG_sequences.fasta -outseq "$EXTRA_DIR"
# If your input FASTA is already proteins, copy instead:
cp all_ARG_sequences.fasta "$EXTRA_DIR"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Activate the Prokka environment
conda activate prokka_env

# Loop through each filtered contig file in the directory
for CONTIG_FILE in "$CONTIGS_DIR"/*_filtered_contigs.fa; do
    # Extract the sample name from the filename (e.g., SSR1_filtered_contigs.fa -> sample1)
    SAMPLE_NAME=$(basename "$CONTIG_FILE" | sed 's/_filtered_contigs.fa//')

    # Define the output directory for this sample
    SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/$SAMPLE_NAME"

    # Create output directory if it doesn't exist
    mkdir -p "$SAMPLE_OUTPUT_DIR"

    # Run Prokka with the --metagenome option and the translated ARG protein database
    prokka --outdir "$SAMPLE_OUTPUT_DIR" \
           --prefix "ARG_$SAMPLE_NAME" \
           --metagenome "$CONTIG_FILE" \
           --force \
           --proteins "$EXTRA_DIR"

    # Check if Prokka ran successfully
    if [ $? -eq 0 ]; then
        echo "Prokka annotation completed for $SAMPLE_NAME"
    else
        echo "Prokka annotation failed for $SAMPLE_NAME"
    fi
done
