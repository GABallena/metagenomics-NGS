\
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Generate per-sample Krona HTML plots from a BIOM-converted TSV abundance table.
#
# Expected input format (typical "biom convert --to-tsv"):
#   - Line 1: optional comment starting with '#'
#   - Header line: sample columns + a final "taxonomy" column
#   - Data lines: first column = feature/OTU ID, sample columns = numeric abundances, last column = taxonomy path
#
# Output:
#   - One Krona HTML per sample

usage() {
  cat <<'EOF'
Usage:
  biom_abundance_to_krona_per_sample.sh [options]

Options:
  -i, --input FILE         Input TSV (default: kraken_output/abundance.tsv)
  -o, --outdir DIR         Output directory for HTML (default: krona_html)
  -t, --tmpdir DIR         Temporary directory (default: krona_tmp)
  -m, --min FLOAT          Minimum abundance to include (default: 0.1)
      --delim STR          Taxonomy delimiter (default: ';')
      --quiet              Suppress ktImportText stdout/stderr
  -h, --help               Show help

Notes:
  - Requires 'ktImportText' in PATH (KronaTools).
  - Abundances can be counts or percentages; Krona accepts numeric weights.
EOF
}

INPUT_TSV="kraken_output/abundance.tsv"
OUTPUT_DIR="krona_html"
TEMP_DIR="krona_tmp"
MIN_ABUNDANCE="0.1"
DELIM=";"
QUIET=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT_TSV="$2"; shift 2;;
    -o|--outdir) OUTPUT_DIR="$2"; shift 2;;
    -t|--tmpdir) TEMP_DIR="$2"; shift 2;;
    -m|--min) MIN_ABUNDANCE="$2"; shift 2;;
    --delim) DELIM="$2"; shift 2;;
    --quiet) QUIET=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown option: $1" >&2; usage; exit 2;;
  esac
done

if [[ ! -f "$INPUT_TSV" ]]; then
  echo "[ERROR] Input TSV not found: $INPUT_TSV" >&2
  exit 1
fi

if ! command -v ktImportText >/dev/null 2>&1; then
  echo "[ERROR] ktImportText not found in PATH. Install KronaTools." >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

# Determine header line: if first line starts with '#', header is line 2; else line 1.
HEADER_LINE=1
first_line="$(head -n 1 "$INPUT_TSV" || true)"
if [[ "${first_line:0:1}" == "#" ]]; then
  HEADER_LINE=2
fi

# Extract sample columns from the header line: columns 2..(NF-1)
mapfile -t samples < <(
  awk -F'\t' -v hl="$HEADER_LINE" '
    NR==hl {
      for (i=2; i<NF; i++) {
        print i "\t" $i
      }
      exit
    }' "$INPUT_TSV"
)

if [[ ${#samples[@]} -eq 0 ]]; then
  echo "[ERROR] Could not parse sample columns from header (line $HEADER_LINE) in: $INPUT_TSV" >&2
  exit 1
fi

echo "[INFO] Parsed ${#samples[@]} sample column(s) from $INPUT_TSV (header line $HEADER_LINE)."
echo "[INFO] Writing per-sample Krona HTML files to: $OUTPUT_DIR/"

for entry in "${samples[@]}"; do
  col_num="$(cut -f1 <<<"$entry")"
  sample_name="$(cut -f2- <<<"$entry")"
  # sanitize sample_name for filenames
  safe_name="$(echo "$sample_name" | tr ' /' '__' | tr -cd '[:alnum:]_.-')"
  input_file="${TEMP_DIR}/${safe_name}.krona.input"
  out_html="${OUTPUT_DIR}/${safe_name}.html"

  # Build Krona input: weight \t level1 \t level2 ...
  awk -F'\t' -v hl="$HEADER_LINE" -v col="$col_num" -v min="$MIN_ABUNDANCE" -v delim="$DELIM" '
    function trim(s) { sub(/^[ \t\r\n]+/, "", s); sub(/[ \t\r\n]+$/, "", s); return s }
    BEGIN { OFS="\t" }
    NR < hl { next }          # skip comment lines above header
    NR == hl { next }         # skip header line
    {
      v = $col + 0
      if (v < min) next
      tax = $(NF)
      tax = trim(tax)
      if (tax == "" || tax == "NA") next

      # split taxonomy on delimiter; allow either ";" or "; " styles
      n = split(tax, a, delim)
      printf "%s", v

      for (i=1; i<=n; i++) {
        lvl = trim(a[i])
        # remove common prefixes like k__, p__, etc
        sub(/^[A-Za-z]__/, "", lvl)
        # drop empty placeholders
        if (lvl == "" || lvl == "NA") continue
        printf "\t%s", lvl
      }
      printf "\n"
    }
  ' "$INPUT_TSV" > "$input_file"

  if [[ ! -s "$input_file" ]]; then
    echo "[WARN] No entries >= ${MIN_ABUNDANCE} for sample '${sample_name}'. Skipping."
    continue
  fi

  if [[ "$QUIET" -eq 1 ]]; then
    ktImportText "$input_file" -o "$out_html" >/dev/null 2>&1
  else
    ktImportText "$input_file" -o "$out_html"
  fi

  echo "[OK] Generated: $out_html"
done

echo "[DONE] Krona HTML outputs in: $OUTPUT_DIR/"
