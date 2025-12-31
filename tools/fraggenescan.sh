#!/bin/bash
set -euo pipefail

# Default parameters (override with FRAGGENESCAN_BIN env var)
FRAGGENESCAN="${FRAGGENESCAN_BIN:-run_FragGeneScan.pl}"
TRAIN_MODEL="complete"  # For complete genes
THREAD_NUM=4

# Verify FragGeneScan exists
if ! command -v "$FRAGGENESCAN" >/dev/null 2>&1 && [[ ! -x "$FRAGGENESCAN" ]]; then
    echo "Error: FragGeneScan binary not found. Set FRAGGENESCAN_BIN or add it to PATH."
    exit 1
fi

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT_FILE="$2"
            shift 2
            ;;
        --output)
            OUTPUT_BASE="$2"
            shift 2
            ;;
        --sample)
            SAMPLE_ID="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

# Validate inputs
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Input file not found: $INPUT_FILE"
    exit 1
fi

if [[ -z "$OUTPUT_BASE" ]]; then
    echo "Output base path not specified"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_BASE")
mkdir -p "$OUTPUT_DIR"

# Run FragGeneScan with explicit path
echo "Running FragGeneScan on sample $SAMPLE_ID"
"$FRAGGENESCAN" \
    -genome="$INPUT_FILE" \
    -out="$OUTPUT_BASE" \
    -complete=1 \
    -train="$TRAIN_MODEL" \
    -thread="$THREAD_NUM" || {
        echo "FragGeneScan failed with exit code $?"
        exit 1
    }

# Check if output files were created
if [[ ! -f "${OUTPUT_BASE}.faa" ]]; then
    echo "FragGeneScan failed to create output file: ${OUTPUT_BASE}.faa"
    exit 1
fi

echo "FragGeneScan completed successfully"
