#!/usr/bin/env bash
# -------------------------------------------------------------------------
# Portfolio-ready note:
# - No hard-coded local paths.
# - Designed to be run from any directory with -i INPUT_DIR provided.
# -------------------------------------------------------------------------

set -euo pipefail

# Rarefy paired-end FASTQs sequentially (no GNU parallel).
# Works with: SAMPLE_R1_trimmed.fq.gz / SAMPLE_R2_trimmed.fq.gz
#
# Usage:
#   ./02_rarefy_reads_seq.sh \
#     -i /path/to/trimmed_reads \
#     -g "1e5,3e5,1e6,3e6,1e7" \
#     -r 5 \
#     -o work/rarefied \
#     -S 1337 \
#     --engine reformat     # or: --engine rasusa
#
# Tips for big inputs (200M pairs):
#   export BBMAP_JAVA_OPTS="-Xmx16g"   # for reformat.sh
#   Keep gzip level low for speed (ZLEVEL=1).

IN=""
GRID=""
REPS=3
OUT="work/rarefied"
SEED=42
ENGINE="reformat"      # or "rasusa"
ZLEVEL="${ZLEVEL:-1}"  # gzip compression level (1=fast)
DRYRUN=0               # if 1, skip heavy work and just plan

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) IN="$2"; shift 2;;
    -g) GRID="$2"; shift 2;;
    -r) REPS="$2"; shift 2;;
    -o) OUT="$2"; shift 2;;
    -S) SEED="$2"; shift 2;;
    --engine) ENGINE="$2"; shift 2;;
    --dry-run) DRYRUN=1; shift 1;;
    *) echo "[ERR] Unknown arg: $1" >&2; exit 1;;
  esac
done

[[ -d "$IN" ]] || { echo "[ERR] -i INPUT_DIR not found: $IN" >&2; exit 1; }
[[ -n "$GRID" ]] || { echo "[ERR] -g depth grid required (e.g., \"1e5,3e5,1e6\")" >&2; exit 1; }
mkdir -p "$OUT"

# Tool checks (skip in dry-run)
if [[ "$ENGINE" == "reformat" ]]; then
  if (( ! DRYRUN )); then
    command -v reformat.sh >/dev/null 2>&1 || { echo "[ERR] reformat.sh not found. Install bbmap." >&2; exit 1; }
  fi
  # Enable pigz if available (no pigt flag to keep compatibility)
  export BBMAP_PIGZ=true
  export BBMAP_GZIP_LEVEL="$ZLEVEL"
elif [[ "$ENGINE" == "rasusa" ]]; then
  if (( ! DRYRUN )); then
    command -v rasusa >/dev/null 2>&1 || { echo "[ERR] rasusa not found. Install rasusa or use --engine reformat." >&2; exit 1; }
  fi
else
  echo "[ERR] --engine must be reformat or rasusa" >&2; exit 1
fi

# Helpers
mate()   { echo "$1" | sed -E 's/_R1_trimmed\.fq\.gz$/_R2_trimmed.fq.gz/'; }
sname()  { basename "$1" | sed -E 's/_R1_trimmed\.fq\.gz$//'; }
to_int() {                    # parse "1e6" or "3,000,000" -> integer (pure bash)
  local s="$1"
  # strip commas and whitespace
  s="${s//,/}"
  s="$(printf '%s' "$s" | tr -d '[:space:]')"
  if [[ $s =~ ^([0-9]+)[eE]([0-9]+)$ ]]; then
    local base="${BASH_REMATCH[1]}"
    local exp="${BASH_REMATCH[2]}"
    local p=1
    local i
    for ((i=0; i<exp; i++)); do p=$((p*10)); done
    echo $(( base * p ))
  else
    # assume integer string
    echo "${s%%.*}"
  fi
}
count_reads() {               # count reads from R1 fastq.gz without awk
  local f="$1"
  if command -v seqkit >/dev/null 2>&1; then
    seqkit stats -Ta "$f" | sed -n '2p' | cut -f4
  else
    # count lines and divide by 4
    local lines
    lines=$(zcat "$f" | wc -l | tr -d ' ')
    echo $(( lines / 4 ))
  fi
}

# Enumerate samples
mapfile -t R1S < <(find "$IN" -maxdepth 1 -name "*_R1_trimmed.fq.gz" \( -type f -o -xtype f \) -print | sort)
[[ ${#R1S[@]} -gt 0 ]] || { echo "[ERR] No *_R1_trimmed.fq.gz files in $IN" >&2; exit 1; }

# Prepare summary
SUMMARY="$OUT/summary.tsv"
echo -e "sample\tdepth\trep\tseed\twritten_pairs" > "$SUMMARY"

# Parse depth grid once
IFS=',' read -ra DEPTHS <<< "$GRID"

# Process samples sequentially
for R1 in "${R1S[@]}"; do
  R2="$(mate "$R1")"
  [[ -s "$R2" ]] || { echo "[ERR] Mate file not found for $R1" >&2; exit 1; }
  SAMPLE="$(sname "$R1")"

  echo "[INFO] Sample: $SAMPLE"
  if (( DRYRUN )); then
    TOTAL=999999999
    echo "[INFO] (dry-run) Skipping count; using TOTAL=$TOTAL"
  else
    TOTAL="$(count_reads "$R1")"
    echo "       Total pairs in R1: $TOTAL"
  fi

  for d in "${DEPTHS[@]}"; do
    TGT="$(to_int "$d")"
    # Cap at total to avoid upsampling errors
    if (( TGT > TOTAL )); then TGT="$TOTAL"; fi
    echo "  - Depth target: $TGT"

    for rep in $(seq 1 "$REPS"); do
      seed=$(( SEED + 131071*rep + 17*TGT ))
      outdir="$OUT/${SAMPLE}/D${TGT}/rep${rep}"
      mkdir -p "$outdir"
      O1="$outdir/${SAMPLE}_R1.fastq.gz"
      O2="$outdir/${SAMPLE}_R2.fastq.gz"
      LOG="$outdir/rarefy.log"

      echo "      rep $rep (seed=$seed) -> $O1"
      if (( DRYRUN )); then
        echo "      [dry-run] $ENGINE would write $O1 and $O2 (target=$TGT)" | tee -a "$LOG"
        WRITTEN="$TGT"
        echo -e "${SAMPLE}\t${TGT}\t${rep}\t${seed}\t${WRITTEN}" >> "$SUMMARY"
        continue
      fi

      if [[ "$ENGINE" == "reformat" ]]; then
        # Exact-by-count, pair-safe
        # Optional JVM mem: export BBMAP_JAVA_OPTS="-Xmx16g"
        reformat.sh in1="$R1" in2="$R2" out1="$O1" out2="$O2" \
          samplereadstarget="$TGT" sampleseed="$seed" tossjunk ow=t \
          ziplevel="$ZLEVEL" pigz=t unpigz=t 2> "$LOG"
      else
        # rasusa exact-by-count, pair-safe
        rasusa -i "$R1" "$R2" -o "$O1" "$O2" -n "$TGT" -s "$seed" > "$LOG" 2>&1
      fi

      # Verify/write summary
      if [[ -s "$O1" ]]; then
        # count written reads without awk
        lines=$(zcat "$O1" | wc -l | tr -d ' ')
        WRITTEN=$(( lines / 4 ))
      else
        WRITTEN="0"
        echo "[WARN] Missing output for $SAMPLE D$TGT rep$rep (see $LOG)" >&2
      fi
      echo -e "${SAMPLE}\t${TGT}\t${rep}\t${seed}\t${WRITTEN}" >> "$SUMMARY"
    done
  done
done

echo "[OK] Rarefaction done."
echo "Head of summary:"
head -n 20 "$SUMMARY" | column -t
