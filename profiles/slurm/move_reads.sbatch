#!/bin/bash
#SBATCH --job-name=stage_reads
#SBATCH --output=logs/stage_reads_%j.out
#SBATCH --error=logs/stage_reads_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal

set -euo pipefail

# ----------------------------
# Usage:
#   sbatch move_reads.sbatch /path/to/source_reads /path/to/project
# ----------------------------

SRC_READS="${1:?ERROR: source read directory required}"
WORKDIR="${2:-$PWD}"

DEST_READS="${WORKDIR}/raw_reads"

mkdir -p "${WORKDIR}/logs"
mkdir -p "${DEST_READS}"

echo "Source reads: ${SRC_READS}"
echo "Destination: ${DEST_READS}"
echo "Started at: $(date)"

rsync -avh --progress "${SRC_READS}/"*.fastq.gz "${DEST_READS}/"

echo "Finished at: $(date)"
