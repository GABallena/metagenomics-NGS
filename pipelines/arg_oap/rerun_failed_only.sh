#!/usr/bin/env bash
# Re-run only samples listed in arg_oap_work/failures.txt
# - Uses conda run to find args_oap even if env isn't activated
# - Hardened TMPDIR + ulimit
# - Decompress+retry if gz streaming is flaky
# - Shard sanity test uses seqtk/reformat.sh/seqkit head (streaming)

set -Eeuo pipefail

BASE="${BASE:-$PWD}"
READS_DIR="${READS_DIR:-$BASE/trimmed_reads}"
WORK_BASE="${WORK_BASE:-$BASE/arg_oap_work}"
FAIL="${FAIL:-$WORK_BASE/failures.txt}"

THREADS="${THREADS:-4}"
ULIMIT_NOFILE="${ULIMIT_NOFILE:-4096}"
ENV_NAME="${ENV_NAME:-args_oap}"

ts(){ date +"%Y-%m-%d %H:%M:%S"; }

# --- pick a runner for args_oap ---
RUN=""
if command -v args_oap >/dev/null 2>&1; then
  RUN=""  # already on PATH
else
  # prefer conda run
  if command -v conda >/dev/null 2>&1; then
    RUN="conda run -n ${ENV_NAME}"
  elif [[ -n "${CONDA_EXE:-}" ]]; then
    RUN="${CONDA_EXE} run -n ${ENV_NAME}"
  else
    echo "ERROR: args_oap not found and conda is unavailable."
    echo "TIP: either 'conda activate ${ENV_NAME}' or install conda, then rerun."
    exit 1
  fi
fi

# --- environment hardening ---
export TMPDIR="${TMPDIR:-$WORK_BASE/tmp}"; mkdir -p "$TMPDIR"
ulimit -n "$ULIMIT_NOFILE" || true

# --- sample list from failures.txt ---
if [[ ! -s "$FAIL" ]]; then
  echo "No $FAIL found or it's empty."; exit 0
fi
mapfile -t SAMPLES < <(awk -F'[ =]+' '/SAMPLE=/{print $4}' "$FAIL" | sort -u)
[[ ${#SAMPLES[@]} -gt 0 ]] || { echo "No samples parsed from $FAIL"; exit 0; }
echo "Will retry: ${SAMPLES[*]}"

# --- helpers ---
ensure_dirs(){
  local S="$1"
  mkdir -p "$WORK_BASE/sample_in/$S" "$WORK_BASE/samples/$S/logs"
  rm -f "$WORK_BASE/samples/$S/.stage1.done" "$WORK_BASE/samples/$S/.stage2.done"
}

logpath(){ echo "$WORK_BASE/samples/$1/logs/$1.log"; }

run_stage1_gz(){
  local S="$1" s_in="$WORK_BASE/sample_in/$S" s_out="$WORK_BASE/samples/$S"
  ln -sf "$READS_DIR/${S}_R1_trimmed.fq.gz" "$s_in/"
  ln -sf "$READS_DIR/${S}_R2_trimmed.fq.gz" "$s_in/"
  echo "[RUN ] Stage 1 (gz) :: $S"
  ${RUN} args_oap stage_one -i "$s_in" -o "$s_out" -t "$THREADS" -f fq \
    --e1 1e-10 --e2 3 --id 45 --qcov 0 |& tee "$(logpath "$S")"
  [[ -s "$s_out/metadata.txt" && -s "$s_out/extracted.fa" ]] && touch "$s_out/.stage1.done"
}

decompress_and_retry(){
  local S="$1" s_in="$WORK_BASE/sample_in/$S" s_out="$WORK_BASE/samples/$S"
  echo "[INFO] Decompressing to .fq and retrying :: $S"
  mkdir -p "$WORK_BASE/unzipped"
  if command -v pigz >/dev/null 2>&1; then ZD="pigz -dc"; else ZD="gzip -dc"; fi
  $ZD "$READS_DIR/${S}_R1_trimmed.fq.gz" > "$WORK_BASE/unzipped/${S}_R1_trimmed.fq"
  $ZD "$READS_DIR/${S}_R2_trimmed.fq.gz" > "$WORK_BASE/unzipped/${S}_R2_trimmed.fq"
  ln -sf "$WORK_BASE/unzipped/${S}_R1_trimmed.fq" "$s_in/"
  ln -sf "$WORK_BASE/unzipped/${S}_R2_trimmed.fq" "$s_in/"
  ${RUN} args_oap stage_one -i "$s_in" -o "$s_out" -t "$THREADS" -f fq \
    --e1 1e-10 --e2 3 --id 45 --qcov 0 |& tee -a "$(logpath "$S")"
  [[ -s "$s_out/metadata.txt" && -s "$s_out/extracted.fa" ]] && touch "$s_out/.stage1.done"
}

shard_test(){
  local S="$1" outdir="$WORK_BASE/samples/${S}_shard"
  echo "[TEST] 2M paired-read shard :: $S"
  mkdir -p "$WORK_BASE/shards" "$outdir/logs"
  local R1="$READS_DIR/${S}_R1_trimmed.fq.gz" R2="$READS_DIR/${S}_R2_trimmed.fq.gz"
  if command -v seqtk >/dev/null 2>&1; then
    seqtk sample -s22 "$R1" 2000000 > "$WORK_BASE/shards/${S}_R1_2M.fq"
    seqtk sample -s22 "$R2" 2000000 > "$WORK_BASE/shards/${S}_R2_2M.fq"
  elif command -v reformat.sh >/dev/null 2>&1; then
    reformat.sh in1="$R1" in2="$R2" out1="$WORK_BASE/shards/${S}_R1_2M.fq" \
      out2="$WORK_BASE/shards/${S}_R2_2M.fq" samplereads=2000000 1>/dev/null
  elif command -v seqkit >/dev/null 2>&1; then
    # streaming, low-RAM (biased head but fine for sanity)
    seqkit head -n 2000000 "$R1" > "$WORK_BASE/shards/${S}_R1_2M.fq"
    seqkit head -n 2000000 "$R2" > "$WORK_BASE/shards/${S}_R2_2M.fq"
  else
    echo "[WARN] No sampler found (seqtk/reformat.sh/seqkit). Skipping shard test."
    return 0
  fi
  gzip -f "$WORK_BASE/shards/${S}_R1_2M.fq" "$WORK_BASE/shards/${S}_R2_2M.fq"
  ${RUN} args_oap stage_one -i "$WORK_BASE/shards" -o "$outdir" -t "$THREADS" -f fq \
    --e1 1e-10 --e2 3 --id 45 --qcov 0 |& tee "$outdir/logs/${S}_shard.log" || true
  if [[ -s "$outdir/metadata.txt" && -s "$outdir/extracted.fa" ]]; then
    echo "[OK  ] Shard passed for $S"
  else
    echo "[WARN] Shard did not pass for $S"
  fi
}

run_stage2(){
  local S="$1" s_out="$WORK_BASE/samples/$S"
  echo "[RUN ] Stage 2 :: $S"
  if ${RUN} args_oap stage_two -i "$s_out" -t "$THREADS" \
       --e 1e-7 --id 80 --qcov 75 --length 25 |& tee -a "$(logpath "$S")"; then
    touch "$s_out/.stage2.done"
  else
    echo "[WARN] Stage 2 returned non-zero for $S (see log)."
  fi
}

# --- main loop ---
for S in "${SAMPLES[@]}"; do
  echo -e "\n== [$S] =="
  rm -rf "$WORK_BASE/samples/$S"
  ensure_dirs "$S"

  run_stage1_gz "$S" || true
  if [[ ! -f "$WORK_BASE/samples/$S/.stage1.done" ]]; then
    decompress_and_retry "$S" || true
  fi

  if [[ ! -f "$WORK_BASE/samples/$S/.stage1.done" ]]; then
    echo "[FAIL] $S :: Stage 1 still failing. Errors:"
    L="$(logpath "$S")"
    [[ -f "$L" ]] && (grep -Ein 'ERROR|WARN|Killed|No space|denied|cannot|fail|oom|out of memory|CRC|truncat|segmentation' "$L" || tail -n 120 "$L") || true
    shard_test "$S"
    continue
  fi

  run_stage2 "$S"
done

echo -e "\nDone. Logs in $WORK_BASE/samples/<S>/logs/"
