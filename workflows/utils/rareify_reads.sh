#!/usr/bin/env bash
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

# rareify_reads.sh — Step 10: build rarefaction grid locally (no Slurm)
# Requires: bash, gzip/pigz, and BBTools (reformat.sh). Optional: GNU coreutils, md5sum.
# Usage:
#   chmod +x rareify_reads.sh
#   ./rareify_reads.sh --read-dir trimmed_reads --out-dir rarefied_reads --depths "250000 500000 1000000 2000000 4000000" --threads 8
#
# Notes:
# - Depths are in **read pairs** (e.g., 1,000,000 means 1M pairs). reformat.sh expects reads (not pairs),
#   so we multiply by 2 internally.
# - Pairing is preserved. Seeds are deterministic per (sample, depth).
# - A manifest CSV is produced at: OUT_DIR/rarefaction_manifest.csv

set -euo pipefail

READ_DIR="trimmed_reads"
OUT_DIR="rarefied_reads"
DEPTHS="250000 500000 1000000 2000000 4000000 8000000"
THREADS=4
DRYRUN=0

# --- Parse args ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --read-dir) READ_DIR="$2"; shift 2;;
    --out-dir) OUT_DIR="$2"; shift 2;;
    --depths) DEPTHS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --dry-run) DRYRUN=1; shift 1;;
    --help|-h)
      sed -n '1,60p' "$0"; exit 0;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

# --- Checks ---
if ! command -v reformat.sh >/dev/null 2>&1; then
  echo "[ERROR] BBTools reformat.sh not found. Install via: conda install -c bioconda bbmap" >&2
  exit 2
fi
mkdir -p "$OUT_DIR"
MANIFEST="${OUT_DIR}/rarefaction_manifest.csv"

# Write header if manifest doesn't exist
if [[ ! -f "$MANIFEST" ]]; then
  echo "sample,depth_pairs,out_R1,out_R2,seed,reads_available_pairs" > "$MANIFEST"
fi

# --- Helper: compute available pairs from R1 fastq ---
get_pairs() {
  local fq1="$1"
  # Fast count via gzip -cd and wc -l; robust but reads whole file.
  # If pigz exists, prefer it for speed.
  if command -v pigz >/dev/null 2>&1; then
    pigz -cd "$fq1" | awk 'END{print NR/4}'
  else
    gzip -cd "$fq1" | awk 'END{print NR/4}'
  fi
}

# --- Helper: derive deterministic seed from sample+depth (first 8 hex of md5 -> int) ---
mk_seed() {
  local s="$1"
  local hex=$(printf "%s" "$s" | md5sum | cut -c1-8)
  # bash base-16 to base-10
  echo $((16#$hex))
}

# --- Discover R1 files ---
shopt -s nullglob
mapfile -t R1S < <(ls -1 "${READ_DIR}"/*_R1_*fastq.gz "${READ_DIR}"/*_R1_*fq.gz "${READ_DIR}"/*_R1.fastq.gz "${READ_DIR}"/*_R1.fq.gz 2>/dev/null || true)

if [[ ${#R1S[@]} -eq 0 ]]; then
  echo "[ERROR] No R1 FASTQs found in ${READ_DIR}. Expected patterns like *_R1_*.fq.gz" >&2
  exit 3
fi

echo "[INFO] Found ${#R1S[@]} samples in ${READ_DIR}"
echo "[INFO] Depth grid (pairs): ${DEPTHS}"

for R1 in "${R1S[@]}"; do
  # Infer R2 by replacing _R1 with _R2 in common patterns
  base="$(basename "$R1")"
  dir="$(dirname "$R1")"

  if [[ -f "${R1/_R1_/_R2_}" ]]; then
    R2="${R1/_R1_/_R2_}"
    sample="${base%%_R1_*}"
  elif [[ -f "${R1/_R1./_R2.}" ]]; then
    R2="${R1/_R1./_R2.}"
    sample="${base%%_R1.*}"
  else
    echo "[WARN] Could not find R2 for $R1; skipping."
    continue
  fi

  # Compute available pairs from R1
  avail_pairs=$(get_pairs "$R1")
  echo "[INFO] Sample ${sample}: available pairs = ${avail_pairs}"

  # Per-sample output directory
  s_out="${OUT_DIR}/${sample}"
  mkdir -p "$s_out"

  for depth in $DEPTHS; do
    if (( depth > avail_pairs )); then
      echo "[WARN] ${sample}: requested depth ${depth} > available ${avail_pairs}; skipping this depth."
      continue
    fi

    seed=$(mk_seed "${sample}_${depth}")
    out1="${s_out}/${sample}.d${depth}_R1.fq.gz"
    out2="${s_out}/${sample}.d${depth}_R2.fq.gz"

    echo "[RUN] ${sample} @ ${depth} pairs (seed=${seed}) -> ${out1}, ${out2}"
    echo "${sample},${depth},${out1},${out2},${seed},${avail_pairs}" >> "$MANIFEST"

    if (( DRYRUN == 1 )); then
      continue
    fi

    # reformat.sh expects number of reads, not pairs -> multiply by 2
    reads_target=$(( depth * 2 ))

    # Execute sampling with preserved pairing and deterministic seed
    reformat.sh \
      in="$R1" in2="$R2" \
      out1="$out1" out2="$out2" \
      samplereadstarget="$reads_target" sampleseed="$seed" \
      ow=t t="$THREADS" gzip=6 1>"${out1}.log" 2>&1

    # Basic sanity: verify outputs exist
    if [[ ! -s "$out1" || ! -s "$out2" ]]; then
      echo "[ERROR] Missing output for ${sample} depth ${depth}" >&2
      exit 4
    fi
  done
done

echo "[DONE] Rarefaction complete. Manifest: ${MANIFEST}"
