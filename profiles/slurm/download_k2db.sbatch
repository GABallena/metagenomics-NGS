#!/bin/bash
#SBATCH --job-name=download_kraken2_db
#SBATCH --output=logs/k2db_%j.out
#SBATCH --error=logs/k2db_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal

set -euo pipefail

# ----------------------------
# Usage:
#   sbatch download_k2db.sbatch /path/to/project
# ----------------------------

WORKDIR="${1:-$PWD}"
DB_DIR="${WORKDIR}/Kraken"

mkdir -p "${WORKDIR}/logs"
mkdir -p "${DB_DIR}"
cd "${DB_DIR}"

# Example public Kraken2 DB (users should choose appropriate version)
KRAKEN_DB_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_2023.tar.gz"
DB_TAR="$(basename "${KRAKEN_DB_URL}")"

echo "Downloading Kraken2 DB to: ${DB_DIR}"
echo "Started at: $(date)"

wget -c "${KRAKEN_DB_URL}"

echo "Extracting database..."
tar -xzf "${DB_TAR}"

echo "Finished at: $(date)"
echo "Contents:"
ls -lh
