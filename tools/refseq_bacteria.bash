#!/bin/bash

# Set directories for database storage and output
DB_DIR="databases/refseq_bacteria_protein"
FASTA_OUTPUT="databases/refseq_bacteria_protein/refseq_bacteria.faa"

# Create directory if it doesn't exist
mkdir -p $DB_DIR

# Step 1: Dynamically find the number of available files by listing them from the FTP server
BASE_URL="https://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/"
echo "Fetching the list of available files from the FTP server..."
available_files=$(curl -s $BASE_URL | grep -oP 'bacteria\.wp_protein\.\d+\.protein\.faa\.gz' | sort -u)

if [ -z "$available_files" ]; then
    echo "No files found to download. Exiting."
    exit 1
fi

# Step 2: Download the RefSeq Bacteria protein sequences using aria2c
echo "Downloading the RefSeq Bacteria protein sequences with aria2c..."
for file in $available_files; do
    URL="${BASE_URL}${file}"
    
    # Check if the file already exists
    if [ -f "$DB_DIR/$file" ]; then
        echo "$file already exists. Skipping download."
    else
        aria2c -x 16 -s 16 -d $DB_DIR $URL  # Download with parallel connections
        if [ $? -ne 0 ]; then
            echo "Download of $file failed. Exiting."
            exit 1
        fi
    fi
done

# Step 3: Verify the download (optional but recommended)
echo "Verifying downloaded files..."
for file in $DB_DIR/*.faa.gz; do
    if [[ -f $file ]]; then
        echo "$file downloaded successfully."
    else
        echo "Error: $file did not download correctly."
        exit 1
    fi
done

# Step 4: Extract the .faa.gz files
echo "Extracting the RefSeq Bacteria protein sequences..."
for gz_file in $DB_DIR/*.faa.gz; do
    gunzip $gz_file
done

# Step 5: Combine all extracted FASTA files into a single FASTA file
echo "Combining all FASTA files into a single file..."
cat $DB_DIR/*.faa > $FASTA_OUTPUT

echo "FASTA file generated at $FASTA_OUTPUT."
echo "RefSeq Bacteria protein database download and processing complete."
