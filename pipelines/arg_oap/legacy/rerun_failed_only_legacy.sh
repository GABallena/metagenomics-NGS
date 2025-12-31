# rerun_failed_only.sh
set -euo pipefail
BASE="$PWD"
FAIL="$BASE/arg_oap_work/failures.txt"
SAMPLES=$(awk -F'[ =]+' '/SAMPLE=/{print $3}' "$FAIL" | sort -u)

echo "Will retry: $SAMPLES"
for S in $SAMPLES; do
  rm -rf "$BASE/arg_oap_work/samples/$S"
  mkdir -p "$BASE/arg_oap_work/samples/$S/logs"
done

# hardened env
export TMPDIR="$BASE/arg_oap_work/tmp"; mkdir -p "$TMPDIR"
ulimit -n 4096 || true

# lower threads a bit to reduce peak load
THREADS=4 RESUME_ONLY=1 READS_DIR="$BASE/trimmed_reads" WORK_BASE="$BASE/arg_oap_work" \
  ./run_arg_oap_local.sh
