\
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Generate Krona HTML plots from Kraken2 .k2report files at multiple ranks.
# Produces separate single-rank Krona plots (not a full hierarchical Krona).

usage() {
  cat <<'EOF'
Usage:
  kraken_report_to_krona_by_rank.sh [options]

Options:
  -i, --indir DIR         Directory with *.k2report files (default: kraken_output)
  -o, --outdir DIR        Base output directory (default: krona_by_rank)
  -t, --tmpdir DIR        Temporary directory (default: krona_by_rank_input)
  -m, --min FLOAT         Minimum percent threshold (column 1) (default: 0.1)
      --counts FIELD      Which count field to use: clade|taxon (default: clade)
      --quiet             Suppress ktImportText stdout/stderr
  -h, --help              Show help

Output layout:
  <outdir>/<rank>/<sample>.html
EOF
}

KRAKEN_DIR="kraken_output"
BASE_OUT="krona_by_rank"
TMP_BASE="krona_by_rank_input"
MIN_PCT="0.1"
COUNTS_FIELD="clade"
QUIET=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--indir) KRAKEN_DIR="$2"; shift 2;;
    -o|--outdir) BASE_OUT="$2"; shift 2;;
    -t|--tmpdir) TMP_BASE="$2"; shift 2;;
    -m|--min) MIN_PCT="$2"; shift 2;;
    --counts) COUNTS_FIELD="$2"; shift 2;;
    --quiet) QUIET=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown option: $1" >&2; usage; exit 2;;
  esac
done

if ! command -v ktImportText >/dev/null 2>&1; then
  echo "[ERROR] ktImportText not found in PATH. Install KronaTools." >&2
  exit 1
fi

count_col=2
if [[ "$COUNTS_FIELD" == "taxon" ]]; then
  count_col=3
elif [[ "$COUNTS_FIELD" != "clade" ]]; then
  echo "[ERROR] Invalid --counts value: $COUNTS_FIELD (use clade|taxon)" >&2
  exit 2
fi

declare -A RANKS=(
  [phylum]="P"
  [class]="C"
  [order]="O"
  [family]="F"
  [genus]="G"
)

files=( "$KRAKEN_DIR"/*.k2report )
if [[ ${#files[@]} -eq 0 ]]; then
  echo "[ERROR] No .k2report files found in: $KRAKEN_DIR" >&2
  exit 1
fi

mkdir -p "$BASE_OUT" "$TMP_BASE"

echo "[INFO] Generating per-rank Krona HTML for ${#files[@]} report(s) (min%=${MIN_PCT}, counts=${COUNTS_FIELD})."

for rank_name in "${!RANKS[@]}"; do
  rank_code="${RANKS[$rank_name]}"
  outdir="${BASE_OUT}/${rank_name}"
  tmpdir="${TMP_BASE}/${rank_name}"
  mkdir -p "$outdir" "$tmpdir"

  for file in "${files[@]}"; do
    sample="$(basename "$file" .k2report)"
    safe_sample="$(echo "$sample" | tr ' /' '__' | tr -cd '[:alnum:]_.-')"
    input_file="${tmpdir}/${safe_sample}.krona.input"
    out_html="${outdir}/${safe_sample}.html"

    awk -F'\t' -v min="$MIN_PCT" -v rank="$rank_code" -v cc="$count_col" '
      BEGIN { OFS="\t" }
      $1+0 >= min && $4 == rank {
        name = $6
        gsub(/^[ ]+/, "", name)
        if (name == "" || name == "NA") next
        print $(cc), name
      }
    ' "$file" > "$input_file"

    if [[ ! -s "$input_file" ]]; then
      continue
    fi

    if [[ "$QUIET" -eq 1 ]]; then
      ktImportText "$input_file" -o "$out_html" >/dev/null 2>&1
    else
      ktImportText "$input_file" -o "$out_html"
    fi
  done

  echo "[OK] Rank '${rank_name}' -> ${outdir}/"
done

echo "[DONE] Krona HTML outputs in: $BASE_OUT/"
