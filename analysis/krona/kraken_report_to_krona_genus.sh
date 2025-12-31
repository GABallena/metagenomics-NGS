\
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Generate Krona HTML plots from Kraken2 .k2report files, filtered to a single taxonomic rank (default: genus).
#
# Kraken2 report columns (tab-delimited):
#   1: % of reads in clade
#   2: reads in clade
#   3: reads assigned directly to this taxon
#   4: rank code (e.g., P,C,O,F,G,S)
#   5: NCBI taxid
#   6+: scientific name (often with leading spaces)
#
# Output:
#   - One Krona HTML per sample

usage() {
  cat <<'EOF'
Usage:
  kraken_report_to_krona_genus.sh [options]

Options:
  -i, --indir DIR         Directory with *.k2report files (default: kraken_output)
  -o, --outdir DIR        Output directory for HTML (default: krona_plots)
  -t, --tmpdir DIR        Temporary directory (default: krona_input)
  -m, --min FLOAT         Minimum percent threshold (column 1) (default: 0.1)
  -r, --rank CODE         Rank code to include (default: G)
      --counts FIELD      Which count field to use: clade|taxon (default: clade)
      --quiet             Suppress ktImportText stdout/stderr
  -h, --help              Show help

Notes:
  - Requires 'ktImportText' in PATH (KronaTools).
EOF
}

KRAKEN_DIR="kraken_output"
KRONA_HTML_DIR="krona_plots"
KRONA_INPUT_DIR="krona_input"
MIN_PCT="0.1"
RANK_CODE="G"
COUNTS_FIELD="clade"  # clade -> column 2, taxon -> column 3
QUIET=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--indir) KRAKEN_DIR="$2"; shift 2;;
    -o|--outdir) KRONA_HTML_DIR="$2"; shift 2;;
    -t|--tmpdir) KRONA_INPUT_DIR="$2"; shift 2;;
    -m|--min) MIN_PCT="$2"; shift 2;;
    -r|--rank) RANK_CODE="$2"; shift 2;;
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

mkdir -p "$KRONA_INPUT_DIR" "$KRONA_HTML_DIR"

count_col=2
if [[ "$COUNTS_FIELD" == "taxon" ]]; then
  count_col=3
elif [[ "$COUNTS_FIELD" != "clade" ]]; then
  echo "[ERROR] Invalid --counts value: $COUNTS_FIELD (use clade|taxon)" >&2
  exit 2
fi

files=( "$KRAKEN_DIR"/*.k2report )
if [[ ${#files[@]} -eq 0 ]]; then
  echo "[ERROR] No .k2report files found in: $KRAKEN_DIR" >&2
  exit 1
fi

echo "[INFO] Converting ${#files[@]} Kraken2 report(s) to Krona HTML (rank=${RANK_CODE}, min%=${MIN_PCT}, counts=${COUNTS_FIELD})."

for file in "${files[@]}"; do
  sample="$(basename "$file" .k2report)"
  safe_sample="$(echo "$sample" | tr ' /' '__' | tr -cd '[:alnum:]_.-')"
  input_file="${KRONA_INPUT_DIR}/${safe_sample}.krona.input"
  out_html="${KRONA_HTML_DIR}/${safe_sample}.html"

  awk -F'\t' -v min="$MIN_PCT" -v rank="$RANK_CODE" -v cc="$count_col" '
    BEGIN { OFS="\t" }
    $1+0 >= min && $4 == rank {
      name = $6
      gsub(/^[ ]+/, "", name)
      if (name == "" || name == "NA") next
      print $(cc), name
    }
  ' "$file" > "$input_file"

  if [[ ! -s "$input_file" ]]; then
    echo "[WARN] No entries matched for sample '${sample}' (rank=${RANK_CODE}, min%=${MIN_PCT}). Skipping."
    continue
  fi

  if [[ "$QUIET" -eq 1 ]]; then
    ktImportText "$input_file" -o "$out_html" >/dev/null 2>&1
  else
    ktImportText "$input_file" -o "$out_html"
  fi

  echo "[OK] Generated: $out_html"
done

echo "[DONE] Krona HTML outputs in: $KRONA_HTML_DIR/"
