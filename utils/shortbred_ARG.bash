#!/bin/bash

# Define directories and file paths
TARGET_DIR="merged_reads"  # Directory containing input folders
COMBINED_DB="databases/ARG+UniProt.fasta"  # Path to the combined database file
ARG_PROTEINS="databases/arg_proteins.fasta"  # Path to the ARG proteins file
UNIPROT_DB="databases/uniprot_sprot.fasta"  # Path to the UniProt proteins file
FILTERED_ARG_DIR="filtered_ARGs"  # Directory to store filtered ARG proteins
MARKERS_DIR="markers_output"  # Directory to store ShortBRED marker files
FILTERED_ARG_PROTEINS="$FILTERED_ARG_DIR/filtered_arg_proteins.fasta"  # Output file for filtered ARG proteins
FILTER_KEYWORDS_FILE="filter_keywords.txt"  # File containing filter keywords
HIGHCONF_ARG_PROTEINS="databases/highconfidence_ARGs.fasta"  # Output file for high-confidence ARG proteins
TSV_FILE="header_mapping.tsv"  # TSV file to store original headers and their corresponding numbers
CLEANED_READS="trimmed_reads"  # Trimmed or cleaned reads
CONTIGS_DIR="assembled_contigs"  # Directory to store assembled contigs

# Create necessary directories
mkdir -p "$FILTERED_ARG_DIR" || { echo "Failed to create $FILTERED_ARG_DIR"; exit 1; }
mkdir -p "$MARKERS_DIR" || { echo "Failed to create $MARKERS_DIR"; exit 1; }
mkdir -p "$TARGET_DIR" || { echo "Failed to create $TARGET_DIR"; exit 1; }
mkdir -p "$CONTIGS_DIR" || { echo "Failed to create $CONTIGS_DIR"; exit 1; }

# Move all files with '_paired_fastq.gz' in their name to the target directory
echo "Copying trimmed/cleaned reads to $TARGET_DIR"
find "$CLEANED_READS" -name '*_paired.fastq.gz' -exec cp {} "$TARGET_DIR"/ \; || { echo "Failed to copy files to $TARGET_DIR"; exit 1; }

# Unzip all .fastq.gz files in the target directory in parallel
echo "Unzipping files in parallel..."
find "$TARGET_DIR" -name "*.fastq.gz" | xargs -P 12 gunzip || { echo "Failed to unzip files"; exit 1; }

# Assemble reads into contigs using MEGAHIT
echo "Running MEGAHIT assembly..."
for fastq_file1 in "$TARGET_DIR"/*_1.fastq; do
    fastq_file2="${fastq_file1/_1.fastq/_2.fastq}"
    base_name=$(basename "$fastq_file1" _1.fastq)
    output_contigs="$CONTIGS_DIR/${base_name}_filtered_contigs.fa"
    
    ./run_megahit.sh "$fastq_file1" "$fastq_file2" "$output_contigs" || { echo "MEGAHIT assembly failed for $base_name"; exit 1; }
    
    echo "MEGAHIT assembly completed for $base_name"
done

# Initialize TSV file
echo -e "Original_Header\tNew_Header" > "$TSV_FILE"
counter=60000000001

# Function to simplify headers
simplify_headers() {
    while read -r line; do
        if [[ $line == ">"* ]]; then
            original_header=$(echo "$line" | cut -d '|' -f2)
            new_header=">$counter"
            echo -e "$original_header\t$new_header" >> "$TSV_FILE"
            echo "$new_header"
            counter=$((counter + 1))
        else
            echo "$line"
        fi
    done
}

# Step 1: Filter the arg_proteins.fasta file using grep and filter keywords
echo "Filtering ARGs based on criteria..."
grep -B1 -i -f "$FILTER_KEYWORDS_FILE" "$ARG_PROTEINS" | grep -v "^--$" > "$FILTERED_ARG_PROTEINS.tmp" || { echo "Filtering ARGs failed"; exit 1; }

# Step 2: Exclude filtered ARGs to create a high-confidence ARGs file
awk 'NR==FNR{if($1 ~ /^>/) a[substr($1,2)]; next} /^>/{f=!(substr($1,2) in a)} f' "$FILTERED_ARG_PROTEINS.tmp" "$ARG_PROTEINS" | simplify_headers > "$HIGHCONF_ARG_PROTEINS" || { echo "Failed to create high-confidence ARG proteins file"; exit 1; }

# Clean up the temporary file
rm "$FILTERED_ARG_PROTEINS.tmp"

# Simplify headers for the combined DB
cat "$HIGHCONF_ARG_PROTEINS" "$UNIPROT_DB" | simplify_headers > "$COMBINED_DB" || { echo "Failed to combine and simplify headers"; exit 1; }

# Activate ShortBRED environment
conda activate shortbred_env || { echo "Failed to activate shortbred_env"; exit 1; }

# Build the ShortBRED database from the high-confidence ARG file
shortbred_identify.py --goi "$HIGHCONF_ARG_PROTEINS" --ref "$COMBINED_DB" --markers "$MARKERS_DIR/highconfidence_ARG_markers.fasta" --clustid 0.95 --threads 12 || { echo "Failed to build ShortBRED markers"; exit 1; }

# Run ShortBRED Quantify for high-confidence ARGs
for contigs_file in "$CONTIGS_DIR"/*.fa; do
    base_name=$(basename "$contigs_file" .fa)
    shortbred_quantify.py --markers "$MARKERS_DIR/highconfidence_ARG_markers.fasta" --wgs "$contigs_file" --results "$MARKERS_DIR/${base_name}_highconfidence_ARG_abundances.txt" --threads 12 || { echo "ShortBRED quantify failed for $contigs_file"; exit 1; }
done

# Print completion message
echo "ShortBRED analysis completed."
