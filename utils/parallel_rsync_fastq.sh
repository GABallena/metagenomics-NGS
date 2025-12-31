#!/bin/bash
set -euo pipefail

# --------------------------------------------
# Parallel rsync utility for paired-end FASTQ files
#
# Usage:
#   ./parallel_rsync_fastq.sh /path/to/local_reads user@remote:/path/to/destination
# --------------------------------------------

SRC_DIR="${1:?ERROR: source directory required}"
DEST="${2:?ERROR: destination (user@host:/path) required}"

LOG_DIR="$(pwd)/upload_logs"
TMP_DIR="$(pwd)/temp_chunks"
GROUPS=4

# Clean up old artifacts
rm -rf "${TMP_DIR}" "${LOG_DIR}"
mkdir -p "${TMP_DIR}" "${LOG_DIR}"
cd "${TMP_DIR}"

echo "Source directory: ${SRC_DIR}"
echo "Destination: ${DEST}"
echo "Parallel groups: ${GROUPS}"
echo "Started at: $(date)"

# STEP 1: Extract unique sample IDs (strip _1.fq.gz / _2.fq.gz)
ls -1 "${SRC_DIR}"/*.fq.gz \
  | sed -E 's/_([12])\.fq\.gz$//' \
  | sort -u > sample_ids.txt

# STEP 2: Split samples evenly into groups
TOTAL=$(wc -l < sample_ids.txt)
PER_GROUP=$(( (TOTAL + GROUPS - 1) / GROUPS ))
split -l "${PER_GROUP}" sample_ids.txt chunk_

# STEP 3: Parallel upload per group
i=0
for CHUNK in chunk_*; do
    ((i++))
    LOG_FILE="${LOG_DIR}/upload_group_${i}.log"

    (
        while IFS= read -r SAMPLE; do
            echo "Uploading sample: ${SAMPLE}"
            rsync -avz --partial --append-verify \
                "${SRC_DIR}/${SAMPLE}_1.fq.gz" "${DEST}/"
            rsync -avz --partial --append-verify \
                "${SRC_DIR}/${SAMPLE}_2.fq.gz" "${DEST}/"
        done < "${CHUNK}"
    ) > "${LOG_FILE}" 2>&1 &
done

wait

echo "All uploads completed"

# STEP 4: Summary
for LOG in "${LOG_DIR}"/*.log; do
    if grep -qi "rsync error" "${LOG}"; then
        echo "Upload failed: $(basename "${LOG}")"
    else
        echo "Upload successful: $(basename "${LOG}")"
    fi
done

echo "Finished at: $(date)"
