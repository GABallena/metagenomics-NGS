#!/usr/bin/env bash
#
# Portfolio-safe copy: paths/identifiers generalized; databases not included.
#
# 00_harmonize_argdb_local.sh
# Sequence-centric harmonization of ARG databases (CARD, MEGARes, ResFinder,
# AMRFinder-NT if present, and SARG if ARG_db/sarg-curation exists)

set -euo pipefail

DBROOT="${1:-}"
if [[ -z "$DBROOT" || ! -d "$DBROOT" ]]; then
  echo "[ERR] Provide the ARG_db root dir (contains card_db, megares3_db, resfinder_db, amrfinder_db, sarg-curation)." >&2
  exit 1
fi

# Auto-detect SARG inside ARG_db
SARGDIR=""
for d in "sarg-curation" "SARG-curation" "sarg_curation"; do
  [[ -d "$DBROOT/$d" ]] && SARGDIR="$DBROOT/$d" && break
done

MIN_LEN="${MIN_LEN:-90}"
WRAP_COLS=60
OUT="work/argdb_harmonized"; TMP="$OUT/tmp"
mkdir -p "$OUT" "$TMP"

need(){ command -v "$1" >/dev/null 2>&1 || { echo "[ERR] Missing $1 (install seqkit, cd-hit)."; exit 1; }; }
need seqkit; need cd-hit-est; need sha256sum; need awk; need grep; need sed; need stat

echo "[INFO] ARG_db:   $DBROOT"
[[ -n "$SARGDIR" ]] && echo "[INFO] SARG dir: $SARGDIR"

POOL="$TMP/pool.fa"; : > "$POOL"

gate_and_append() {
  local src="$1"; local file="$2"
  [[ -s "$file" ]] || return 0
  # quick protein guard: if <60% ACGTN, skip
  local sample raw frac
  sample=$(head -n 200 "$file" | tr -d '\n[:space:]' | tr -cd 'ACGTNacgtn' | wc -c || true)
  raw=$(head -n 200 "$file" | tr -d '\n[:space:]' | tr -cd 'A-Za-z' | wc -c || true)
  if [[ "${raw:-0}" -gt 0 ]]; then
    frac=$(( 100 * sample / raw ))
    if [[ "$frac" -lt 60 ]]; then
      echo "[SKIP] $src: $(basename "$file") looks non-DNA (≤60% ACGTN)."
      return 0
    fi
  fi
  seqkit seq -t DNA -u "$file" \
  | awk -v SRC="$src" -v F="$(basename "$file")" -v MINLEN="$MIN_LEN" -v WRAP="$WRAP_COLS" '
      BEGIN{RS=">"; ORS=""; OFS=" "}
      NR==1{next}
      {
        split($0, L, "\n"); hdr=L[1]; seq=""
        for(i=2;i<=length(L);i++) seq=seq L[i]
        gsub(/[^ACGTN]/,"",seq)
        if(length(seq) < MINLEN) next
        split(hdr,a," "); acc=a[1]
        printf(">%s SRC=%s FILE=%s RAW=%s\n", acc, SRC, F, hdr)
        for(i=1;i<=length(seq);i+=WRAP) print substr(seq,i,WRAP) "\n"
      }' >> "$POOL"
  echo "[OK] + ${src}: $(basename "$file") (DNA, len≥$MIN_LEN)"
}

# --- Collect nucleotide FASTAs ---
CARD="$DBROOT/card_db"
gate_and_append "CARD" "$CARD/nucleotide_fasta_protein_homolog_model.fasta"
gate_and_append "CARD" "$CARD/nucleotide_fasta_protein_variant_model.fasta"
gate_and_append "CARD" "$CARD/nucleotide_fasta_protein_knockout_model.fasta"
gate_and_append "CARD" "$CARD/nucleotide_fasta_protein_overexpression_model.fasta"
gate_and_append "CARD" "$CARD/nucleotide_fasta_rRNA_gene_variant_model.fasta"

MEG="$DBROOT/megares3_db"
gate_and_append "MEGARES" "$MEG/megares_database_v3.00.fasta"

RES="$DBROOT/resfinder_db"
gate_and_append "RESFINDER" "$RES/all.fsa"

AMR="$DBROOT/amrfinder_db"
if [[ -d "$AMR" ]]; then
  # include only nucleotide files that pass the DNA gate
  while IFS= read -r -d '' f; do
    gate_and_append "AMRFINDER" "$f"
  done < <(find "$AMR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) -print0)
fi

if [[ -n "$SARGDIR" ]]; then
  for sargf in "$SARGDIR"/sarg*.fa; do
    [[ -s "$sargf" ]] && gate_and_append "SARG" "$sargf" || true
  done
fi

grep -q '^>' "$POOL" || { echo "[ERR] No DNA sequences passed the gate."; exit 1; }

# --- Input manifest ---
{
  echo -e "file\tbytes\tsha256"
  find "$DBROOT" -maxdepth 2 -type f \( -name "*.fa" -o -name "*.fsa" -o -name "*.fasta" \) \
    | sort | while read -r f; do
      [[ -s "$f" ]] || continue
      sz=$(stat -c%s "$f"); sh=$(sha256sum "$f" | awk '{print $1}')
      echo -e "${f}\t${sz}\t${sh}"
    done
  if [[ -n "$SARGDIR" ]]; then
    for f in "$SARGDIR"/sarg*.fa; do
      [[ -s "$f" ]] || continue
      sz=$(stat -c%s "$f"); sh=$(sha256sum "$f" | awk '{print $1}')
      echo -e "${f}\t${sz}\t${sh}"
    done
  fi
} > "$OUT/downloads_manifest.tsv"

# --- Normalize to sequence-hash IDs + synonyms ---
STD_FA="$TMP/pool.std.fa"
SYN="$OUT/synonyms.tsv"
echo -e "hash\tlength\tsource_db\tsource_file\taccession\traw_header" > "$SYN"

awk -v SYN="$SYN" -v WRAP="$WRAP_COLS" -v MINLEN="$MIN_LEN" '
  BEGIN{RS=">"; OFS="\t"; ORS=""}
  NR==1{next}
  {
    split($0, L, "\n"); hdr=L[1]; seq=""
    for(i=2;i<=length(L);i++) seq=seq L[i]
    gsub(/[^ACGTN]/,"",seq)
    if(length(seq) < MINLEN) next
    match(hdr,/SRC=([^ ]+)/,m1); src=(m1[1]!=""?m1[1]:"NA")
    match(hdr,/FILE=([^ ]+)/,m2); fil=(m2[1]!=""?m2[1]:"NA")
    raw=hdr; sub(/.*RAW=/,"",raw)
    acc=hdr; sub(/ .*/,"",acc)
    cmd="printf \"%s\" \"" seq "\" | sha256sum | awk \047{print $1}\047"
    cmd | getline h; close(cmd)
    print h"\t"length(seq)"\t"src"\t"fil"\t"acc"\t"raw"\n" >> SYN
    printf(">%s\n", h)
    for(i=1;i<=length(seq);i+=WRAP) printf("%s\n", substr(seq,i,WRAP))
  }' "$POOL" > "$STD_FA"

# --- CD-HIT 100% ---
NR_FA="$OUT/nonredundant.fa"
CLSTR="$OUT/cdhit_100.clstr"
cd-hit-est -i "$STD_FA" -o "$NR_FA" -c 1.0 -aS 1.0 -aL 1.0 -G 0 -g 1 -d 0 -T "$(nproc)" -M 64000
mv "${NR_FA}.clstr" "$CLSTR"

CLMAP="$OUT/cluster_map.tsv"
echo -e "orig_hash\trep_hash" > "$CLMAP"
awk '
  /^>Cluster/ {rep=""}
  /\*$/       {match($0,/>[^ ]+/,m); rep=substr(m[0],2)}
  /at /       {match($0,/>[^ ]+/,m); h=substr(m[0],2); if(rep=="") rep=h; print h"\t"rep}
' "$CLSTR" >> "$CLMAP"

# --- Manifests & summary ---
{
  echo -e "file\tbytes\tsha256"
  for f in "$NR_FA" "$SYN" "$CLMAP" "$CLSTR" "$STD_FA" "$POOL"; do
    [[ -s "$f" ]] || continue
    sz=$(stat -c%s "$f"); sh=$(sha256sum "$f" | awk '{print $1}')
    echo -e "${f}\t${sz}\t${sh}"
  done
} > "$OUT/outputs_manifest.tsv"

echo "[DONE] Nonredundant FASTA: $NR_FA (records: $(grep -c '^>' "$NR_FA"))"
echo "[DONE] Synonym map:        $SYN     (rows:    $(($(wc -l < "$SYN")-1)))"
echo "[DONE] Cluster map:        $CLMAP   (rows:    $(($(wc -l < "$CLMAP")-1)))"
echo "[DONE] CD-HIT clusters:    $CLSTR"

