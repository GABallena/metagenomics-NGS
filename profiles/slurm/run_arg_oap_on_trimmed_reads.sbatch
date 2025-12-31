#!/usr/bin/env bash
#SBATCH --job-name=arg_oap_batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# run_arg_oap_on_trimmed_reads.sbatch
# Portfolio-safe SBATCH wrapper for running ARG-OAP over many paired-end samples.
#
# Assumptions (edit as needed):
#   - Paired reads live in: trimmed_reads/
#   - File naming: <SAMPLE>_R1_*.fq.gz and <SAMPLE>_R2_*.fq.gz
#   - The actual run logic is implemented in:
#       pipelines/arg_oap/run_arg_oap_local.sh
#
# Typical usage:
#   sbatch profiles/slurm/run_arg_oap_on_trimmed_reads.sbatch
#
# Optional overrides:
#   READS_DIR=/path/to/trimmed_reads WORK_BASE=/path/to/work sbatch ...
#
set -Eeuo pipefail

mkdir -p logs

BASE="${BASE:-$PWD}"
export READS_DIR="${READS_DIR:-$BASE/trimmed_reads}"
export WORK_BASE="${WORK_BASE:-$BASE/arg_oap_work}"

# If your cluster uses a shared scratch, point TMPDIR there:
export TMPDIR="${TMPDIR:-$BASE/tmp}"
mkdir -p "$TMPDIR"

# Threading: default to SLURM allocation
export THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"

echo "[INFO] BASE=$BASE"
echo "[INFO] READS_DIR=$READS_DIR"
echo "[INFO] WORK_BASE=$WORK_BASE"
echo "[INFO] THREADS=$THREADS"
echo "[INFO] Host=$(hostname)  Date=$(date)"

bash "$BASE/pipelines/arg_oap/run_arg_oap_local.sh"
