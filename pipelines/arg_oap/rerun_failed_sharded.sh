#!/usr/bin/env bash
# Rerun only the failed samples (from arg_oap_work/failures.txt) using sharding.
# For each sample S:
#   1) seqkit split2 by read pairs (low RAM)
#   2) Stage 1 per shard (modest threads)
#   3) Merge Stage 1 outputs
#   4) Single Stage 2 on merged outputs
# Produces .stage1.done / .stage2.done in arg_oap_work/samples/<S>

set -Eeuo pipefail

# ---------------- user knobs ----------------
BASE="${BASE:-$PWD}"
READS_DIR="${READS_DIR:-$BASE/trimmed_reads}"
WORK_BASE="${WORK_BASE:-$BASE/arg_oap_work}"
FAIL_FILE="${FAIL_FILE:-$WORK_BASE/failures.txt}"

# Sharding & threads
PAIRS_PER_SHARD="${PAIRS_PER_SHARD:-15000000}"  # ~15M pairs/shard; tune 10–20M
THREADS_STAGE1="${THREADS_STAGE1:-4}"           # per-shard threads
THREADS_STAGE2="${THREADS_STAGE2:-8}"           # final Stage 2 threads

# ARG-OAP params (same as your wrapper)
E1="${E1:-1e-10}"; E2="${E2:-3}"; ID1="${ID1:-45}"; QCOV1="${QCOV1:-0}"
EVAL2="${EVAL2:-1e-7}"; PID2="${PID2:-80}"; QCOV2="${QCOV2:-75}"; ALNLEN2="${ALNLEN2:-25}"

# Env/runner
ENV_NAME="${ENV_NAME:-args_oap}"
ULIMIT_NOFILE="${ULIMIT_NOFILE:-4096}"
# -------------------------------------------

die(){ echo "ERROR: $*" >&2; exit 1; }
ts(){ date +"%Y-%m-%d %H:%M:%S"; }

need(){ command -v "$1" >/dev/null 2>&1 || die "Missing dependency: $1"; }

# runner for args_oap (works even if env not activated)
pick_runner(){
  if command -v args_oap >/dev/null 2>&1; then
    RUN=""
  elif command -v conda >/dev/null 2>&1; then
    RUN="conda run -n ${ENV_NAME}"
  elif [[ -n "${CONDA_EXE:-}" ]]; then
    RUN="${CONDA_EXE} run -n ${ENV_NAME}"
  else
    die "args_oap not found and conda unavailable. Install/activate and retry."
  fi
  echo "$RUN"
}

main() {
  [[ -s "$FAIL_FILE" ]] || die "No failures.txt at $FAIL_FILE (or empty)."
  mapfile -t SAMPLES < <(awk -F'[ =]+' '/SAMPLE=/{print $4}' "$FAIL_FILE" | sort -u)
  [[ ${#SAMPLES[@]} -gt 0 ]] || die "No SAMPLE= entries parsed from $FAIL_FILE."

  need seqkit
  [[ -d "$READS_DIR" ]] || die "Reads dir not found: $READS_DIR"
  mkdir -p "$WORK_BASE/tmp" "$WORK_BASE/shards" "$WORK_BASE/sample_in" "$WORK_BASE/samples"
  export TMPDIR="${TMPDIR:-$WORK_BASE/tmp}"
  ulimit -n "$ULIMIT_NOFILE" || true

  RUN="$(pick_runner)"
  ZD="gzip -dc"; command -v pigz >/dev/null 2>&1 && ZD="pigz -dc"

  echo "[$(ts)] Retrying (sharded): ${SAMPLES[*]}"
  echo "  PAIRS_PER_SHARD=$PAIRS_PER_SHARD  THREADS_STAGE1=$THREADS_STAGE1  THREADS_STAGE2=$THREADS_STAGE2"
  echo

  for S in "${SAMPLES[@]}"; do
    echo "== [$S] =="
    R1="$READS_DIR/${S}_R1_trimmed.fq.gz"
    R2="$READS_DIR/${S}_R2_trimmed.fq.gz"
    [[ -f "$R1" && -f "$R2" ]] || { echo "[SKIP] Missing reads for $S"; continue; }

    # Clean any previous attempt
    rm -rf "$WORK_BASE/samples/$S" "$WORK_BASE/samples/${S}.part_"* "$WORK_BASE/shards/$S"
    mkdir -p "$WORK_BASE/samples/$S/logs"

    SHARD_DIR="$WORK_BASE/shards/$S"
    mkdir -p "$SHARD_DIR"

    echo "[$(ts)]  Sharding to ~$PAIRS_PER_SHARD pairs/part..."
    seqkit split2 -1 "$R1" -2 "$R2" -s "$PAIRS_PER_SHARD" -O "$SHARD_DIR" -e .gz --force

    mapfile -t PARTS < <(find "$SHARD_DIR" -type f -name "*_R1_*part_*.fq.gz" | sort)
    if [[ ${#PARTS[@]} -eq 0 ]]; then
      echo "[FAIL] No shard parts produced for $S — skipping."
      continue
    fi

    # Stage 1 per shard
    ok_shards=0
    for R1P in "${PARTS[@]}"; do
      R2P="${R1P/_R1_/_R2_}"
      [[ -f "$R2P" ]] || die "Missing paired shard for $R1P"

      TAG=$(basename "$R1P" | sed 's/.*\(part_[0-9][0-9][0-9]\).*/\1/')
      S_IN_PART="$WORK_BASE/sample_in/${S}.${TAG}"
      S_OUT_PART="$WORK_BASE/samples/${S}.${TAG}"
      mkdir -p "$S_IN_PART" "$S_OUT_PART/logs"

      ln -sf "$R1P" "$S_IN_PART/"
      ln -sf "$R2P" "$S_IN_PART/"

      echo "[$(ts)]   - Stage 1 :: $S $TAG"
      set +e
      ${RUN} args_oap stage_one \
        -i "$S_IN_PART" -o "$S_OUT_PART" -t "$THREADS_STAGE1" -f fq \
        --e1 "$E1" --e2 "$E2" --id "$ID1" --qcov "$QCOV1" >> "$S_OUT_PART/logs/${S}.${TAG}.log" 2>&1
      rc=$?
      set -e
      if [[ $rc -eq 0 && -s "$S_OUT_PART/metadata.txt" && -s "$S_OUT_PART/extracted.fa" ]]; then
        ok_shards=$((ok_shards+1))
      else
        echo "         [WARN] shard $TAG failed (see $S_OUT_PART/logs/${S}.${TAG}.log)"
      fi
    done

    if [[ $ok_shards -eq 0 ]]; then
      echo "[FAIL] All shards failed for $S — investigate per-shard logs."
      continue
    fi

    # Merge Stage 1 outputs into samples/<S>
    MERGED_DIR="$WORK_BASE/samples/$S"
    : > "$MERGED_DIR/metadata.txt"
    : > "$MERGED_DIR/extracted.fa"

    first=1
    for OUTP in $(find "$WORK_BASE/samples" -maxdepth 1 -type d -name "${S}.part_*" | sort); do
      META="$OUTP/metadata.txt"; EXT="$OUTP/extracted.fa"
      [[ -s "$META" && -s "$EXT" ]] || continue
      if [[ $first -eq 1 ]]; then
        cat "$META" >> "$MERGED_DIR/metadata.txt"; first=0
      else
        awk 'FNR==1 && NR!=1 {next} {print}' "$META" >> "$MERGED_DIR/metadata.txt"
      fi
      cat "$EXT" >> "$MERGED_DIR/extracted.fa"
    done

    if [[ ! -s "$MERGED_DIR/extracted.fa" ]]; then
      echo "[FAIL] Merged extracted.fa is empty for $S — cannot proceed to Stage 2."
      continue
    fi
    touch "$MERGED_DIR/.stage1.done"

    echo "[$(ts)]  Stage 2 (merged) :: $S"
    set +e
    ${RUN} args_oap stage_two \
      -i "$MERGED_DIR" -t "$THREADS_STAGE2" \
      --e "$EVAL2" --id "$PID2" --qcov "$QCOV2" --length "$ALNLEN2" \
      >> "$MERGED_DIR/logs/${S}.merged.log" 2>&1
    rc2=$?
    set -e

    if [[ $rc2 -eq 0 ]]; then
      touch "$MERGED_DIR/.stage2.done"
      echo "[OK ] $S complete (merged)."
    else
      echo "[WARN] Stage 2 returned non-zero for $S (see $MERGED_DIR/logs/${S}.merged.log)"
    fi
  done

  echo
  echo "[$(ts)] Done. Check: $WORK_BASE/samples/<S>/.stage2.done"
}

main "$@"
