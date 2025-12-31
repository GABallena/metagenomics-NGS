#!/bin/bash
set -euo pipefail


# Set directories for database storage and output
DB_DIR="databases/refseq_bacteria"
FASTA_OUTPUT="databases/refseq_bacteria/refseq_bacteria.fasta"

# Create directory if it doesn't exist
mkdir -p $DB_DIR

# Step 1: Download the RefSeq Bacteria protein database using aria2 with no file allocation and overwrite enabled
echo "Downloading the RefSeq Bacteria protein database..."
for i in {1..10}; do
    FILE="bacteria.${i}.protein.faa.gz"
    URL="ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/$FILE"
    aria2c --file-allocation=none --allow-overwrite=true -x 16 -s 16 -d $DB_DIR $URL
    if [ ! -f "$DB_DIR/$FILE" ]; then
        echo "Finished downloading all existing segments at bacteria.${i}.protein.faa.gz"
        break
    fi
done

# Step 2: Verify the download (optional but recommended)
echo "Verifying downloaded files..."
for file in $DB_DIR/*.faa.gz; do
    echo "Checking $file..."
    if [[ -f $file ]]; then
        echo "$file downloaded successfully."
    else
        echo "Error: $file did not download correctly."
        exit 1
    fi
done

# Step 3: Extract the .faa.gz files
echo "Extracting the RefSeq Bacteria protein database files..."
for gz_file in $DB_DIR/*.faa.gz; do
    gunzip $gz_file
done

# Step 4: Combine all extracted FASTA files into a single FASTA file
echo "Combining all FASTA files into a single file..."
cat $DB_DIR/*.faa > $FASTA_OUTPUT

echo "FASTA file generated at $FASTA_OUTPUT."

echo "RefSeq Bacteria database download and processing complete."
