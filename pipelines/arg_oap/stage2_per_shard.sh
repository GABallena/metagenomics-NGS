#!/usr/bin/env bash
# Run ARG-OAP Stage 2 per shard directory: samples/<S>.part_XXX
# Usage:
#   WORK_BASE=... THREADS=4 bash stage2_per_shard.sh SAMPLE-22_RUN7 SAMPLE-23_RUN8 ...

set -Eeuo pipefail

WORK_BASE="${WORK_BASE:-$PWD/arg_oap_work}"
THREADS="${THREADS:-4}"
ENV_NAME="${ENV_NAME:-args_oap}"

# pick runner even if env isn't activated
if command -v args_oap >/dev/null 2>&1; then
  RUN=""
elif command -v conda >/dev/null 2>&1; then
  RUN="conda run -n ${ENV_NAME}"
elif [[ -n "${CONDA_EXE:-}" ]]; then
  RUN="${CONDA_EXE} run -n ${ENV_NAME}"
else
  echo "ERROR: args_oap not found and conda unavailable." >&2; exit 1
fi

EVAL2="${EVAL2:-1e-7}"
PID2="${PID2:-80}"
QCOV2="${QCOV2:-75}"
ALNLEN2="${ALNLEN2:-25}"

ulimit -n 4096 || true
export TMPDIR="${TMPDIR:-$WORK_BASE/tmp}"; mkdir -p "$TMPDIR"

ts(){ date +"%Y-%m-%d %H:%M:%S"; }

[[ $# -ge 1 ]] || { echo "Usage: $0 <SAMPLE_ID> [<SAMPLE_ID> ...]"; exit 1; }

for S in "$@"; do
  echo "== [$S] Stage 2 per shard =="
  # find shard outputs from Stage 1
  mapfile -t SHARDS < <(find "$WORK_BASE/samples" -maxdepth 1 -type d -name "${S}.part_*" | sort)
  if [[ ${#SHARDS[@]} -eq 0 ]]; then
    echo "[WARN] No shards found for $S under $WORK_BASE/samples/${S}.part_*"; continue
  fi

  for D in "${SHARDS[@]}"; do
    [[ -d "$D" ]] || continue
    LOG="$D/logs/$(basename "$D").stage2.log"
    MARK="$D/.stage2.done"
    mkdir -p "$(dirname "$LOG")"

    if [[ -f "$MARK" ]]; then
      echo "  [SKIP] $(basename "$D") stage2 already done."
      continue
    fi
    echo "  [RUN ] $(basename "$D") :: $(ts)"
    set +e
    ${RUN} args_oap stage_two \
      -i "$D" -t "$THREADS" \
      --e "$EVAL2" --id "$PID2" --qcov "$QCOV2" --length "$ALNLEN2" \
      >> "$LOG" 2>&1
    rc=$?
    set -e
    if [[ $rc -eq 0 ]]; then
      touch "$MARK"
      echo "  [OK  ] $(basename "$D")"
    else
      echo "  [FAIL] $(basename "$D") (see $LOG)"
    fi
  done
done
