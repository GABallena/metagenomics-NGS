#!/usr/bin/env bash
# rerun_failed_only_v2.sh
# - Re-runs only samples listed in arg_oap_work/failures.txt
# - Hardens env (TMPDIR, ulimit, lower threads)
# - If Stage 1 fails: auto-decompress to .fq and retry
# - If still fails: run a 2M-read shard sanity test
# - Prints key log lines; keeps failures listed in failures.txt

set -Eeuo pipefail
BASE="${BASE:-$PWD}"
FAIL="$BASE/arg_oap_work/failures.txt"
READS_DIR="${READS_DIR:-$BASE/trimmed_reads}"
WORK_BASE="${WORK_BASE:-$BASE/arg_oap_work}"
THREADS="${THREADS:-4}"                 # keep modest
ULIMIT_NOFILE="${ULIMIT_NOFILE:-4096}"

ts(){ date +"%Y-%m-%d %H:%M:%S"; }

# Parse SAMPLE IDs correctly from failures.txt lines like:
# "2025-08-14 ...  SAMPLE=SAMPLE-22_RUN7  STAGE=stage1  LOG=/path"
mapfile -t SAMPLES < <(awk -F'[ =]+' '/SAMPLE=/{print $4}' "$FAIL" | sort -u)
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "No samples parsed from $FAIL"; exit 0
fi
echo "Will retry: ${SAMPLES[*]}"

export TMPDIR="$WORK_BASE/tmp"; mkdir -p "$TMPDIR"
ulimit -n "$ULIMIT_NOFILE" || true

run_stage1(){
  local S="$1"
  local s_in="$WORK_BASE/sample_in/$S"
  local s_out="$WORK_BASE/samples/$S"
  mkdir -p "$s_in" "$s_out/logs"
  rm -f "$s_out/.stage1.done" "$s_out/.stage2.done"

  # point to gz FASTQs
  ln -sf "$READS_DIR/${S}_R1_trimmed.fq.gz" "$s_in/"
  ln -sf "$READS_DIR/${S}_R2_trimmed.fq.gz" "$s_in/"

  echo "[RUN ] Stage 1 (gz) :: $S"
  args_oap stage_one -i "$s_in" -o "$s_out" -t "$THREADS" -f fq \
    --e1 1e-10 --e2 3 --id 45 --qcov 0 \
    |& tee "$s_out/logs/${S}.log"

  # success criteria
  [[ -s "$s_out/metadata.txt" && -s "$s_out/extracted.fa" ]] && { touch "$s_out/.stage1.done"; return 0; }
  return 1
}

retry_decompress_and_stage1(){
  local S="$1"
  local s_in="$WORK_BASE/sample_in/$S"
  local s_out="$WORK_BASE/samples/$S"
  echo "[INFO] Decompressing to .fq and retrying :: $S"
  mkdir -p "$WORK_BASE/unzipped"
  pigz -dc "$READS_DIR/${S}_R1_trimmed.fq.gz" > "$WORK_BASE/unzipped/${S}_R1_trimmed.fq"
  pigz -dc "$READS_DIR/${S}_R2_trimmed.fq.gz" > "$WORK_BASE/unzipped/${S}_R2_trimmed.fq"
  ln -sf "$WORK_BASE/unzipped/${S}_R1_trimmed.fq" "$s_in/"
  ln -sf "$WORK_BASE/unzipped/${S}_R2_trimmed.fq" "$s_in/"

  args_oap stage_one -i "$s_in" -o "$s_out" -t "$THREADS" -f fq \
    --e1 1e-10 --e2 3 --id 45 --qcov 0 \
    |& tee -a "$s_out/logs/${S}.log"

  [[ -s "$s_out/metadata.txt" && -s "$s_out/extracted.fa" ]] && { touch "$s_out/.stage1.done"; return 0; }
  return 1
}

shard_test(){
  local S="$1"
  local s_out="$WORK_BASE/samples/${S}_shard"
  echo "[TEST] 2M paired-read shard :: $S"
  mkdir -p "$WORK_BASE/shards"
  seqkit sample -s 22 -n 2000000 "$READS_DIR/${S}_R1_trimmed.fq.gz" > "$WORK_BASE/shards/${S}_R1_2M.fq"
  seqkit sample -s 22 -n 2000000 "$READS_DIR/${S}_R2_trimmed.fq.gz" > "$WORK_BASE/shards/${S}_R2_2M.fq"
  gzip -f "$WORK_BASE/shards/${S}_R1_2M.fq" "$WORK_BASE/shards/${S}_R2_2M.fq"
  args_oap stage_one -i "$WORK_BASE/shards" -o "$s_out" -t "$THREADS" -f fq \
    --e1 1e-10 --e2 3 --id 45 --qcov 0 \
    |& tee "$s_out/logs/${S}_shard.log" || true
  [[ -s "$s_out/metadata.txt" && -s "$s_out/extracted.fa" ]] && echo "[OK  ] Shard passed for $S" || echo "[WARN] Shard did not pass for $S"
}

for S in "${SAMPLES[@]}"; do
  echo -e "\n== [$S] =="
  # clean outputs for a fresh attempt
  rm -rf "$WORK_BASE/samples/$S"
  mkdir -p "$WORK_BASE/samples/$S/logs" "$WORK_BASE/sample_in/$S"

  if run_stage1 "$S"; then
    echo "[OK  ] $S :: Stage 1 (gz) succeeded"
  elif retry_decompress_and_stage1 "$S"; then
    echo "[OK  ] $S :: Stage 1 (decompressed) succeeded"
  else
    echo "[FAIL] $S :: Stage 1 still failing. Extracting log tail + errors:"
    L="$WORK_BASE/samples/$S/logs/${S}.log"
    grep -Ein 'ERROR|WARN|Killed|No space|denied|cannot|fail|oom|out of memory|CRC|truncat|segmentation' "$L" || tail -n 120 "$L" || true
    shard_test "$S"
  fi

  # If Stage 1 is done, try Stage 2 right away
  if [[ -f "$WORK_BASE/samples/$S/.stage1.done" ]]; then
    echo "[RUN ] Stage 2 :: $S"
    args_oap stage_two -i "$WORK_BASE/samples/$S" -t "$THREADS" \
      --e 1e-7 --id 80 --qcov 75 --length 25 \
      |& tee -a "$WORK_BASE/samples/$S/logs/${S}.log" || true
    [[ -s "$WORK_BASE/samples/$S/whatever_stage2_outputs_are" ]] && touch "$WORK_BASE/samples/$S/.stage2.done" || true
  fi
done

echo -e "\n[NOTE] Re-run complete. See per-sample logs under $WORK_BASE/samples/<S>/logs/"
