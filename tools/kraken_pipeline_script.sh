#!/bin/bash
set -e  # Stop the script on any error

# Set paths manually (modify these as needed)
TRIMMED_READS_DIR="trimmed_reads"  # Path to trimmed reads
KRAKEN2_DB_DIR="Kraken/k2_db"         # Path to Kraken2 database
KRAKEN2_OUTPUT_DIR="Kraken/kraken_output" # Path to Kraken2 output directory
THREADS_KRAKEN=4                          # Number of threads to use

# Debug: Print paths to ensure they're set correctly
echo "Trimmed reads directory: $TRIMMED_READS_DIR"
echo "Kraken2 database directory: $KRAKEN2_DB_DIR"
echo "Kraken2 output directory: $KRAKEN2_OUTPUT_DIR"
echo "Threads for Kraken2: $THREADS_KRAKEN"

# Validate inputs
if [ ! -d "$TRIMMED_READS_DIR" ] || [ ! -d "$KRAKEN2_DB_DIR" ]; then
    echo "Error: One or more required directories do not exist."
    exit 1
fi

if ! command -v kraken2 &> /dev/null; then
    echo "Error: Kraken2 is not available in the current environment."
    exit 1
fi

if ! [[ "$THREADS_KRAKEN" =~ ^[0-9]+$ ]] || [ "$THREADS_KRAKEN" -le 0 ]; then
    echo "Error: THREADS_KRAKEN must be a positive integer."
    exit 1
fi

# Ensure the output directory exists
mkdir -p "$KRAKEN2_OUTPUT_DIR"

# Activate Kraken2 conda environment
echo "Activating Kraken2 conda environment..."
if [[ -n "${CONDA_BASE:-}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
    source "${HOME}/anaconda3/etc/profile.d/conda.sh"
else
    echo "Error: conda initialization script not found. Set CONDA_BASE to your installation."
    exit 1
fi
conda activate "${KRAKEN2_CONDA_ENV:-kraken_env}"

# Run Kraken2 for each sample
echo "Running Kraken2 on trimmed reads..."
for R1_FILE in "$TRIMMED_READS_DIR"/*_1_paired.fastq; do
    base=$(basename "$R1_FILE" _1_paired.fastq)
    R2_FILE="$TRIMMED_READS_DIR/${base}_2_paired.fastq"

    echo "Debug: Processing sample $base"
    echo "  R1 file: $R1_FILE"
    echo "  R2 file: $R2_FILE"
    echo "  Output report: $KRAKEN2_OUTPUT_DIR/${base}.k2report"
    echo "  Output file: $KRAKEN2_OUTPUT_DIR/${base}.kraken2"

    # Ensure both R1 and R2 files exist
    if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
        echo "Warning: Missing paired files for $base. Skipping."
        continue
    fi

    # Run Kraken2
    kraken2 --db "$KRAKEN2_DB_DIR" \
    --threads "$THREADS_KRAKEN" \
    --report "$KRAKEN2_OUTPUT_DIR/${base}.k2report" \
    --paired "$R1_FILE" "$R2_FILE" \
    > "$KRAKEN2_OUTPUT_DIR/${base}.kraken2" 2>> "$KRAKEN2_OUTPUT_DIR/kraken2_errors.log"

    echo "Debug: Kraken2 processing complete for $base"
done


echo "$(date '+%Y-%m-%d %H:%M:%S') - Kraken2 pipeline complete!" >> "$KRAKEN2_OUTPUT_DIR/kraken2_log.txt"
