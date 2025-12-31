#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Download core_nt database using aria2c
# Ensure aria2c is installed: sudo apt-get install aria2

OUTDIR="${OUTDIR:-core_nt_db}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"
SPLIT="${SPLIT:-16}"

# Set the base URL
BASE_URL="https://ftp.ncbi.nlm.nih.gov/blast/db/core_nt."

# Prepare a list of all the URLs (00 to 68)
echo "Generating core_nt download list..."
for i in $(seq -w 0 68); do
    echo "${BASE_URL}${i}.tar.gz"
    echo "${BASE_URL}${i}.tar.gz.md5"
done > core_nt_download_list.txt

# Download using aria2c with 4 concurrent downloads and 16 connections per file
echo "Starting downloads with aria2c..."
aria2c -i core_nt_download_list.txt --max-concurrent-downloads=${MAX_CONCURRENT} --split=${SPLIT} --continue=true

# Verify MD5 checksums
echo "Verifying checksums..."
for i in $(seq -w 0 68); do
    md5sum -c core_nt.${i}.tar.gz.md5
done

echo "All downloads and verifications completed!"
