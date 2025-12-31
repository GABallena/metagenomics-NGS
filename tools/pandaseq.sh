#!/bin/bash
set -euo pipefail


set -eo pipefail

# Check if required commands are available
for cmd in pandaseq seqkit; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Error: Required command '$cmd' not found" >&2
        exit 1
    fi
done

# Parameters for PANDAseq
MIN_OVERLAP="${MIN_OVERLAP:-10}"
QUALITY_THRESHOLD="${QUALITY_THRESHOLD:-0.9}"
THREADS="${THREADS:-1}"

# Verify pandaseq version
PANDASEQ_VERSION=$(pandaseq -v 2>&1 | head -n1)
echo "Using PANDAseq version: $PANDASEQ_VERSION" >&2

# Function for cleanup on error
cleanup() {
    local exit_code=$?
    echo "Cleaning up temporary files..." >&2
    rm -f "$MERGED_READS.tmp" "${READ1%.gz}.tmp" "${READ2%.gz}.tmp"
    exit $exit_code
}

trap cleanup ERR EXIT

# Validate inputs
for var in READ1 READ2 MERGED_READS QUALITY_METRICS LOG; do
    if [ -z "${!var}" ]; then
        echo "Error: Required variable $var is not set" >&2
        exit 1
    fi
done

for file in "$READ1" "$READ2"; do
    if [ ! -f "$file" ]; then
        echo "Error: Input file $file does not exist" >&2
        exit 1
    fi
    if [ ! -s "$file" ]; then
        echo "Error: Input file $file is empty" >&2
        exit 1
    fi
done

# Validate numeric parameters
if ! [[ "$MIN_OVERLAP" =~ ^[0-9]+$ ]] || ! [[ "$QUALITY_THRESHOLD" =~ ^[0-9.]+$ ]]; then
    echo "Error: MIN_OVERLAP and QUALITY_THRESHOLD must be numeric values" >&2
    exit 1
fi

# Create output directories with error checking
for dir in "$(dirname "$MERGED_READS")" "$(dirname "$QUALITY_METRICS")" "$(dirname "$LOG")"; do
    if ! mkdir -p "$dir"; then
        echo "Error: Failed to create directory: $dir" >&2
        exit 1
    fi
done

# Check available memory
AVAILABLE_MEM=$(free -g | awk '/^Mem:/{print $2}')
if [ "$AVAILABLE_MEM" -lt 4 ]; then
    echo "Warning: Less than 4GB RAM available" >&2
fi

# Run PANDAseq with automatic handling of gzipped files
echo "Running PANDAseq..." >&2
if [[ "$READ1" == *.gz ]] && [[ "$READ2" == *.gz ]]; then
    echo "Decompressing gzipped input files..." >&2
    gzip -dc "$READ1" > "${READ1%.gz}.tmp"
    gzip -dc "$READ2" > "${READ2%.gz}.tmp"
    pandaseq -f "${READ1%.gz}.tmp" -r "${READ2%.gz}.tmp" \
             -o "$MIN_OVERLAP" -t "$QUALITY_THRESHOLD" -T "$THREADS" -N -F \
             -w "$MERGED_READS.tmp" 2> "$LOG" || {
        echo "Error: PANDAseq failed. Check $LOG for details" >&2
        exit 1
    }
    rm -f "${READ1%.gz}.tmp" "${READ2%.gz}.tmp"
else
    pandaseq -f "$READ1" -r "$READ2" -o "$MIN_OVERLAP" \
             -t "$QUALITY_THRESHOLD" -T "$THREADS" -N -F \
             -w "$MERGED_READS.tmp" 2> "$LOG" || {
        echo "Error: PANDAseq failed. Check $LOG for details" >&2
        exit 1
    }
fi

# Validate PANDAseq output
if [ ! -s "$MERGED_READS.tmp" ]; then
    echo "Error: PANDAseq failed to produce output" >&2
    exit 1
fi

# Move temporary file to final location
mv "$MERGED_READS.tmp" "$MERGED_READS"

# Generate quality metrics with progress monitoring
echo "Generating quality metrics..." >&2
if [[ "$MERGED_READS" == *.gz ]]; then
    gzip -dc "$MERGED_READS" | seqkit stats > "$QUALITY_METRICS"
else
    seqkit stats "$MERGED_READS" > "$QUALITY_METRICS"
fi

# Final validation
echo "Validating outputs..." >&2
for file in "$MERGED_READS" "$QUALITY_METRICS"; do
    if [ ! -s "$file" ]; then
        echo "Error: Output file $file is empty" >&2
        exit 1
    fi
done

echo "Processing completed successfully" >&2