#!/usr/bin/env bash
# run_arg_oap_local_resume.sh
# - Skips finished samples using .stage1.done / .stage2.done
# - RESUME_ONLY=1 -> only iterate samples missing .stage2.done
# - CONTINUE_ON_FAIL=1 -> don't exit on failure; append to failures.txt
# - Writes logs per-sample and a consolidated failures.txt

set -Eeuo pipefail

# ==== USER VARS (env-overridable) ============================================
READS_DIR="${READS_DIR:-$PWD/trimmed_reads}"
WORK_BASE="${WORK_BASE:-$PWD/arg_oap_work}"
THREADS="${THREADS:-8}"

STAGE1_CMD="${STAGE1_CMD:-args_oap stage_one}"
STAGE2_CMD="${STAGE2_CMD:-args_oap stage_two}"

# stage_one params
FORMAT_EXT="${FORMAT_EXT:-fq}"     # gz auto-detected
E1="${E1:-1e-10}"
E2="${E2:-3}"
ID1="${ID1:-45}"
QCOV1="${QCOV1:-0}"

# stage_two params
EVAL2="${EVAL2:-1e-7}"
PID2="${PID2:-80}"
QCOV2="${QCOV2:-75}"
ALNLEN2="${ALNLEN2:-25}"

# behavior switches
RESUME_ONLY="${RESUME_ONLY:-0}"         # 1 = only samples missing .stage2.done
CONTINUE_ON_FAIL="${CONTINUE_ON_FAIL:-1}" # 1 = don't exit on failure; log it
ULIMIT_NOFILE="${ULIMIT_NOFILE:-4096}"  # open-file cap
# ============================================================================

ts(){ date +"%Y-%m-%d %H:%M:%S"; }

fail_file_init(){
  mkdir -p "$WORK_BASE"
  FAIL_FILE="$WORK_BASE/failures.txt"
  if [[ ! -f "$FAIL_FILE" ]]; then : > "$FAIL_FILE"; fi
  {
    echo "===== $(ts) run begin ====="
  } >> "$FAIL_FILE"
}

record_failure(){ # $1=sample $2=stage $3=logpath $4=note(optional)
  local s="$1" st="$2" logp="$3" note="${4:-}"
  {
    echo "$(ts)  SAMPLE=${s}  STAGE=${st}  LOG=${logp}  ${note}"
  } >> "$FAIL_FILE"
}

pair_for_r1(){
  local r1="$1" r2=""
  r2="${r1/_R1_/_R2_}"
  [[ -f "$r2" ]] || r2="${r1/_R1/_R2}"
  [[ -f "$r2" ]] || return 1
  printf '%s\n' "$r2"
}

run_sample(){
  local r1="$1" fname base r2 s_in s_out log done1 done2
  fname="$(basename "$r1")"
  base="${fname%%_R1_*}"
  r2="$(pair_for_r1 "$r1")" || { echo "[FAIL] $base (pair)"; record_failure "$base" "pair" "(none)" "R2 not found"; return 1; }

  s_in="${WORK_BASE}/sample_in/${base}"
  s_out="${WORK_BASE}/samples/${base}"
  log="${s_out}/logs/${base}.log"
  done1="${s_out}/.stage1.done"
  done2="${s_out}/.stage2.done"

  mkdir -p "$s_in" "$s_out/logs"

  # Symlink inputs
  ln -sf "$r1" "${s_in}/"
  ln -sf "$r2" "${s_in}/"

  echo "== [$base] =="

  # ---- Stage 1 ----
  if [[ -f "$done1" && -s "${s_out}/metadata.txt" && -s "${s_out}/extracted.fa" ]]; then
    echo "[SKIP] Stage 1 already done."
  else
    echo "[RUN ] Stage 1..."
    {
      echo "== $(ts) :: Stage 1 =="
      echo "CMD: $STAGE1_CMD -i $s_in -o $s_out -t $THREADS -f $FORMAT_EXT --e1 $E1 --e2 $E2 --id $ID1 --qcov $QCOV1"
    } >> "$log"
    set +e
    $STAGE1_CMD -i "$s_in" -o "$s_out" -t "$THREADS" -f "$FORMAT_EXT" \
      --e1 "$E1" --e2 "$E2" --id "$ID1" --qcov "$QCOV1" >> "$log" 2>&1
    rc=$?
    set -e
    if [[ $rc -ne 0 || ! -s "${s_out}/metadata.txt" || ! -s "${s_out}/extracted.fa" ]]; then
      echo "[FAIL] $base (Stage 1). See $log"
      record_failure "$base" "stage1" "$log"
      [[ "$CONTINUE_ON_FAIL" == "1" ]] && return 1 || exit 1
    fi
    touch "$done1"
  fi

  # ---- Stage 2 ----
  if [[ -f "$done2" ]]; then
    echo "[SKIP] Stage 2 already done."
    echo "[OK  ] $base"
    return 0
  else
    echo "[RUN ] Stage 2..."
    {
      echo "== $(ts) :: Stage 2 =="
      echo "CMD: $STAGE2_CMD -i $s_out -t $THREADS --e $EVAL2 --id $PID2 --qcov $QCOV2 --length $ALNLEN2"
    } >> "$log"
    set +e
    $STAGE2_CMD -i "$s_out" -t "$THREADS" \
      --e "$EVAL2" --id "$PID2" --qcov "$QCOV2" --length "$ALNLEN2" >> "$log" 2>&1
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
      echo "[FAIL] $base (Stage 2). See $log"
      record_failure "$base" "stage2" "$log"
      [[ "$CONTINUE_ON_FAIL" == "1" ]] && return 1 || exit 1
    fi
    touch "$done2"
  fi

  echo "[OK  ] $base"
}

preflight(){
  [[ -d "$READS_DIR" ]] || { echo "ERROR: Reads dir not found: $READS_DIR" >&2; exit 1; }
  mkdir -p "$WORK_BASE/sample_in" "$WORK_BASE/samples" "$WORK_BASE/tmp"
  # Use our own tmp to avoid /tmp quota issues
  : "${TMPDIR:=$WORK_BASE/tmp}"
  export TMPDIR
  # Raise file handle limit a bit
  ulimit -n "$ULIMIT_NOFILE" || true
}

collect_samples(){
  mapfile -t R1S < <(find "$READS_DIR" -maxdepth 1 -type f -name "*_R1_*trimmed.fq.gz" | sort)
  if [[ "${RESUME_ONLY}" == "1" ]]; then
    local needed=()
    for r1 in "${R1S[@]}"; do
      local b; b="$(basename "$r1")"; b="${b%%_R1_*}"
      [[ -f "$WORK_BASE/samples/$b/.stage2.done" ]] || needed+=("$r1")
    done
    R1S=("${needed[@]}")
    echo "Resuming: ${#R1S[@]} samples need work (missing .stage2.done)."
  else
    echo "Discovered ${#R1S[@]} R1 files to consider."
  fi
}

main(){
  preflight
  fail_file_init
  collect_samples

  local total=${#R1S[@]} failed=0 processed=0
  echo "Running sequentially at $THREADS threads."
  for r1 in "${R1S[@]}"; do
    if run_sample "$r1"; then
      processed=$((processed+1))
    else
      failed=$((failed+1))
      # continue regardless; failure already recorded
    fi
  done

  echo "----- SUMMARY -----"
  echo "Attempted : $total"
  echo "Succeeded : $processed"
  echo "Failed    : $failed"
  echo "Failures file: $WORK_BASE/failures.txt"
  echo "Outputs in  : $WORK_BASE/samples/"
  {
    echo "===== $(ts) run end: attempted=$total ok=$processed failed=$failed ====="
  } >> "$FAIL_FILE"
}

main "$@"
