#!/bin/bash

# Directories
RAW_DIR="raw_reads"
OUT_DIR="trimmed_reads"
REPORT_DIR="fastqc_output"
LOG_DIR="logs/fastp"

# Make output folders
mkdir -p "$OUT_DIR" "$REPORT_DIR" "$LOG_DIR"

# Loop through all R1 reads
for R1 in ${RAW_DIR}/*_1.fq.gz; do
    # Derive sample name (e.g., sample1 from sample1_1.fq.gz)
    sample=$(basename "$R1" | sed 's/_1\.fq\.gz//')

    R2="${RAW_DIR}/${sample}_2.fq.gz"
    OUT1="${OUT_DIR}/${sample}_R1_paired.fq.gz"
    OUT2="${OUT_DIR}/${sample}_R2_paired.fq.gz"
    HTML="${REPORT_DIR}/${sample}_fastp.html"
    JSON="${REPORT_DIR}/${sample}_fastp.json"
    LOG="${LOG_DIR}/${sample}.log"

    echo "Processing $sample ..."
    fastp -i "$R1" -I "$R2" -o "$OUT1" -O "$OUT2" \
          -w 4 --html "$HTML" --json "$JSON" > "$LOG" 2>&1
done

echo "All done!"
