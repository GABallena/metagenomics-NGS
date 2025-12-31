#!/bin/bash

set -eo pipefail

# Check if SPAdes is available
if ! command -v spades.py >/dev/null 2>&1; then
    echo "Error: SPAdes (spades.py) is not available in PATH" >&2
    exit 1
fi

# Get arguments from Snakemake
READ1="$1"
READ2="$2"
OUTPUT_DIR="$3"
OUTPUT_FILE="$4"
THREADS="$5"
LOG_FILE="$6"
MIN_LEN="$7"

# Function for logging
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

# Clean output directory if it exists
if [ -d "${OUTPUT_DIR}" ]; then
    log "Cleaning existing output directory: ${OUTPUT_DIR}"
    rm -rf "${OUTPUT_DIR}"
fi

# Create fresh output directory
mkdir -p "${OUTPUT_DIR}"

# Check memory
MEM_AVAILABLE=$(free -g | awk '/^Mem:/{print $2}')
MEM_LIMIT=$(( MEM_AVAILABLE * 75 / 100 ))
[ "$MEM_LIMIT" -lt 4 ] && MEM_LIMIT=4  # Set minimum 4GB

# Run SPAdes assembly
log "Starting SPAdes assembly"
SPADES_BIN=$(command -v spades.py)
if ! "${SPADES_BIN}" \
    --meta \
    --only-assembler \
    -1 "${READ1}" \
    -2 "${READ2}" \
    -o "${OUTPUT_DIR}" \
    -t "${THREADS}" \
    -m "${MEM_LIMIT}" \
    -k 21,33,55,77,99,127 \
    --phred-offset 33 \
    2> >(tee -a "${LOG_FILE}"); then
    log "Error: SPAdes assembly failed"
    exit 1
fi

# Verify SPAdes output
CONTIGS_FILE="${OUTPUT_DIR}/contigs.fasta"
if [[ ! -f "${CONTIGS_FILE}" ]]; then
    log "Error: SPAdes output not found: ${CONTIGS_FILE}"
    exit 1
fi

# Filter contigs
log "Filtering contigs with minimum length: ${MIN_LEN}"
if ! awk -v min_len="${MIN_LEN}" \
    'BEGIN {RS=">"} 
     NR>1 {seq=$0; gsub(/\n/,"",seq); 
     if (length(seq) >= min_len) 
        print ">"$0}' \
    "${CONTIGS_FILE}" > "${OUTPUT_FILE}"; then
    log "Error: Contig filtering failed"
    exit 1
fi

# Verify output
if [[ ! -s "${OUTPUT_FILE}" ]]; then
    log "Error: Empty output file"
    exit 1
fi

log "Assembly complete: ${OUTPUT_FILE}"