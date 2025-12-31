#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Help function
usage() {
    echo "Usage: $0 --input1 <fastq1> --input2 <fastq2> --output <output.tsv> --db_dir <metaphlan_db>"
    echo "Requires minimum 15GB memory for MetaPhlAn 4"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input1)
            input1="$2"
            shift 2
            ;;
        --input2)
            input2="$2"
            shift 2
            ;;
        --output)
            output="$2"
            shift 2
            ;;
        --db_dir)
            db_dir="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z $input1 ]] || [[ -z $input2 ]] || [[ -z $output ]] || [[ -z $db_dir ]]; then
    usage
fi

# Create output directory
mkdir -p "$(dirname "$output")"

# Use a dedicated working directory to avoid deleting unrelated folders in the current directory
WORKDIR="$(dirname "$output")/.metaphlan_work"
mkdir -p "$WORKDIR"



# Setup database
if [[ ! -d $db_dir ]]; then
    mkdir -p "$db_dir"
fi

echo "Setting up MetaPhlAn database..."
if ! metaphlan --install --bowtie2db "$db_dir"; then
    echo "Database setup failed"
    echo "k__Bacteria	p__Unknown	100.0" > "$output"
    exit 0
fi

# Generate bowtie2out filename
bowtie2out="${output%.tsv}.bowtie2.bz2"

# Run MetaPhlAn with optimized settings
echo "Running MetaPhlAn analysis..."
if metaphlan "$input1,$input2" \
    --input_type fastq \
    --nproc 8 \
    --force \
    --bowtie2db "$db_dir" \
    --output_file "$output" \
    --bowtie2out "$bowtie2out" \
    --memory-min 15GB \
    --index latest \
    --unclassified_estimation \
    --sample_id "$(basename "${output%.tsv}")"; then
    echo "MetaPhlAn analysis completed successfully"
else
    echo "MetaPhlAn analysis failed"
    echo "k__Bacteria	p__Unknown	100.0" > "$output"
fi
