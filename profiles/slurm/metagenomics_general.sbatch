#!/bin/bash
#SBATCH --job-name=metagenomics_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --qos=normal

set -euo pipefail

# ----------------------------
# Usage:
#   sbatch metagenomics_general.sbatch /path/to/project
# If no argument is given, uses current directory.
# ----------------------------

WORKDIR="${1:-$PWD}"

mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}"

echo "Running pipeline in: ${WORKDIR}"
echo "Started at: $(date)"

# Load conda in a portable way
if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
    source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]; then
    source "${HOME}/anaconda3/etc/profile.d/conda.sh"
else
    echo "ERROR: Conda initialization script not found."
    exit 1
fi

snakemake \
  --snakefile metagenomics_general.smk \
  --use-conda \
  --cores "${SLURM_CPUS_PER_TASK}" \
  --printshellcmds \
  --rerun-incomplete \
  --latency-wait 60 \
  --conda-prefix .snakemake/conda

echo "Finished at: $(date)"
