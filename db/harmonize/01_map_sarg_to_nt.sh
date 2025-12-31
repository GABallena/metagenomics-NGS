#!/usr/bin/env bash
#
# Portfolio-safe copy: paths/identifiers generalized; databases not included.
#
# 01_map_sarg_to_nt.sh  (v4: streaming aggregation, low-memory)
# SARG (protein) → nucleotide-hash bridge via tblastn with HSP aggregation.
# Output: <OUTDIR>/sarg_nt_bridge.tsv (+ tmp, logs, manifest)

set -euo pipefail

# ── CONFIG (override via env) ────────────────────────────────────────────────
THREADS="${THREADS:-16}"

# Pair-level filters (AFTER aggregating all HSPs for a query–target pair)
ID="${ID:-70}"           # min length-weighted AA identity (%)
MIN_AA="${MIN_AA:-50}"   # min total aligned AA length (aa)
QCOV="${QCOV:-60}"       # min query coverage by length (%), capped at 100
SCOV="${SCOV:-50}"       # min subject coverage by nt-union (%)
TOP_MARGIN="${TOP_MARGIN:-10}"   # % bitscore margin around top per query

# External sort tuning
SORT_MEM="${SORT_MEM:-50%}"      # e.g., 8G or 50%
SORT_TMP="${SORT_TMP:-/tmp}"     # temp dir with plenty of space

# Optional reciprocal check (requires mmseqs)
RECIP="${RECIP:-0}"

# ── ARGS ─────────────────────────────────────────────────────────────────────
if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <nonredundant.fa> <outdir> <sarg_protein.fa> [more_sarg.fa ...]" >&2
  exit 1
fi
NT_FA="$1"; shift
OUTDIR="$1"; shift
SARG_PROTS=( "$@" )

mkdir -p "$OUTDIR" "$OUTDIR/tmp" "$OUTDIR/logs"
TMP="$OUTDIR/tmp"

need(){ command -v "$1" >/dev/null 2>&1 || { echo "[ERR] Missing $1 on PATH"; exit 1; }; }
need makeblastdb; need tblastn; need seqkit; need sort; need awk
if [[ "$RECIP" == "1" ]]; then command -v mmseqs >/dev/null 2>&1 || { echo "[ERR] Missing mmseqs for RECIP=1"; exit 1; }; fi
[[ -s "$NT_FA" ]] || { echo "[ERR] Missing nucleotide FASTA: $NT_FA"; exit 1; }
for f in "${SARG_PROTS[@]}"; do [[ -s "$f" ]] || { echo "[ERR] Missing SARG protein FASTA: $f"; exit 1; }; done

# ── 1) Merge + de-dup SARG proteins ─────────────────────────────────────────
SARG_MERGED="$TMP/sarg_proteins.merged.fa"; : > "$SARG_MERGED"
for f in "${SARG_PROTS[@]}"; do
  echo "[INFO] Adding SARG protein file: $f"
  cat "$f" >> "$SARG_MERGED"
done
SARG_AA="$TMP/sarg_proteins.dedup.fa"
seqkit seq -t PROTEIN -u "$SARG_MERGED" | seqkit rmdup -s > "$SARG_AA"
echo "[INFO] SARG queries: $(grep -c '^>' "$SARG_AA") (after dedup)"

# ── 2) Build BLAST NT DB (no -parse_seqids; SHA256 ids are long) ────────────
DBPFX="$OUTDIR/nt_hash_db"
if [[ ! -f "${DBPFX}.nin" && ! -f "${DBPFX}.00.nin" ]]; then
  echo "[INFO] Building BLAST nucleotide DB..."
  makeblastdb -in "$NT_FA" -dbtype nucl -out "$DBPFX" > "$OUTDIR/logs/makeblastdb.log" 2>&1
fi

# ── 3) Run tblastn ──────────────────────────────────────────────────────────
RAW_TSV="$TMP/tblastn_raw.tsv"
echo "[INFO] Running tblastn (threads=$THREADS)..."
tblastn -query "$SARG_AA" -db "$DBPFX" \
  -evalue 1e-10 -seg no -max_target_seqs 200 -num_threads "$THREADS" \
  -outfmt "6 qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
  > "$RAW_TSV"
echo "[INFO] tblastn HSP rows: $(wc -l < "$RAW_TSV")"

# ── 4) External sort by (qacc, sseqid) for streaming aggregation ────────────
# Columns: 1=qacc, 2=sseqid
SORTED="$TMP/tblastn_sorted.tsv"
LC_ALL=C sort -S "$SORT_MEM" -T "$SORT_TMP" -t $'\t' -k1,1 -k2,2 "$RAW_TSV" > "$SORTED"

# ── 5) Stream-aggregate HSPs per (qacc,sseqid); write pair-level rows ───────
# Output columns (pair-level):
# qacc sseqid id_lenw cumAA qcov scov sumBits qlen slen qst qen sst sen
AGG_TSV="$TMP/tblastn_agg.tsv"
awk -v OFS="\t" '
  function flush_pair(   id_lenw,qcov,union_nt,scov,i,cs,ce){
    if(pair_seen==0) return
    # length-weighted identity
    id_lenw = (sum_al>0) ? 100.0*sum_idlen/sum_al : 0
    # qcov capped at 100
    qcov = (qlen>0) ? 100.0*sum_al/qlen : 0; if(qcov>100) qcov=100
    # merge subject intervals → union_nt
    n=interval_n
    if(n>0){
      # sort by start
      for(i=2;i<=n;i++){
        s=starts[i]; e=ends[i]; j=i-1
        while(j>=1 && starts[j]>s){ starts[j+1]=starts[j]; ends[j+1]=ends[j]; j-- }
        starts[j+1]=s; ends[j+1]=e
      }
      union_nt=0; cs=starts[1]; ce=ends[1]
      for(i=2;i<=n;i++){
        if(starts[i] <= ce){ if(ends[i]>ce) ce=ends[i] }
        else { union_nt += (ce-cs+1); cs=starts[i]; ce=ends[i] }
      }
      union_nt += (ce-cs+1)
    } else { union_nt=0 }
    scov = (slen>0) ? 100.0*union_nt/slen : 0
    printf "%s\t%s\t%.2f\t%d\t%.2f\t%.2f\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\n",
           cur_q, cur_s, id_lenw, sum_al, qcov, scov, sum_bits, qlen, slen, first_qst, first_qen, first_sst, first_sen
  }
  BEGIN{
    cur_q=""; cur_s="";
    pair_seen=0
    sum_bits=sum_al=sum_idlen=0
    interval_n=0
  }
  {
    q=$1; s=$2
    if(q!=cur_q || s!=cur_s){
      flush_pair()
      # reset for new pair
      cur_q=q; cur_s=s; pair_seen=0
      sum_bits=sum_al=sum_idlen=0
      interval_n=0
      qlen=$14+0; slen=$15+0
      first_qst=$8+0; first_qen=$9+0; first_sst=$10+0; first_sen=$11+0
    }
    pid=$3+0; al=$4+0; qst=$8+0; qen=$9+0; sst=$10+0; sen=$11+0; bs=$13+0
    if(slen<=0 || qlen<=0) next
    pair_seen=1
    sum_bits += bs
    sum_al   += al
    sum_idlen+= pid*al
    # store subject interval normalized
    a=sst; b=sen; if(a>b){t=a;a=b;b=t}
    interval_n++; starts[interval_n]=a; ends[interval_n]=b
  }
  END{ flush_pair() }
' "$SORTED" > "$AGG_TSV"

echo "[INFO] Pair-level rows (aggregated): $(wc -l < "$AGG_TSV")"

# ── 6) Filter aggregated pairs ───────────────────────────────────────────────
FILTERED="$TMP/pairs_filtered.tsv"
awk -v ID="$ID" -v MINAA="$MIN_AA" -v QCOV="$QCOV" -v SCOV="$SCOV" -v OFS="\t" '
  { pid=$3+0; cumAA=$4+0; qcov=$5+0; scov=$6+0; if(pid>=ID && cumAA>=MINAA && qcov>=QCOV && scov>=SCOV) print $0 }
' "$AGG_TSV" > "$FILTERED"

PASS_N=$(wc -l < "$FILTERED" || echo 0)
echo "[INFO] Passing pairs: $PASS_N"

if [[ "$PASS_N" -eq 0 ]]; then
  OUT_TSV="$OUTDIR/sarg_nt_bridge.tsv"
  echo -e "prot_db\tprot_acc\tprot_len_aa\tprot_name\tnt_hash\tnt_len_bp\taa_identity\taa_aln_len\tqcov_aa\tscov_cds\tbitscore\tevalue\trank_status\tnear_top_hits\tdelta_bitscore_to_next\tframe\tstart_nt\tend_nt\treciprocal_ok" > "$OUT_TSV"
  echo "[NOTE] No aggregated pairs met thresholds. Try easing: ID=60 QCOV=50 SCOV=40 MIN_AA=40"
  exit 0
fi

# ── 7) Rank per query by sum(bitscore) and label status ──────────────────────
RANKED="$TMP/pairs_ranked.tsv"
awk -v MARGIN="$TOP_MARGIN" -v OFS="\t" '
  { q=$1; sb=$7+0; if(!(q in best) || sb>best[q]) best[q]=sb; n[q]++ ; store[q,n[q]]=$0 }
  END{
    for(q in n){
      top=best[q]; cut=top*(1.0 - MARGIN/100.0)
      # count near-top
      c=0; second=0
      for(i=1;i<=n[q];i++){ split(store[q,i],a,"\t"); b=a[7]+0; if(b>=cut) c++; if(b>second && b<top) second=b }
      status = (c==1) ? "single_best" : ( (c<=3) ? "tied_best" : "ambiguous" )
      delta  = (second==0) ? "NA" : sprintf("%.2f",(top-second)/top*100.0)
      for(i=1;i<=n[q];i++){
        split(store[q,i],a,"\t")
        print a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[12],a[13],status,c,delta
      }
    }
  }
' "$FILTERED" > "$RANKED"

# ── 8) Optional reciprocal check (mmseqs tBLASTn-like) ───────────────────────
RECIP_COL="reciprocal_ok"; RECIP_FILE="$TMP/reciprocal.tsv"
if [[ "$RECIP" == "1" ]]; then
  echo "[INFO] Reciprocal (mmseqs search-type 3)..."
  QDB="$TMP/mm_qdb"; SDB="$TMP/mm_sdb"; RES="$TMP/mm_res"
  mmseqs createdb "$SARG_AA" "$QDB" >/dev/null
  mmseqs createdb "$NT_FA" "$SDB" >/dev/null
  mmseqs search "$QDB" "$SDB" "$RES" "$TMP" --threads "$THREADS" --search-type 3 >/dev/null
  mmseqs convertalis "$QDB" "$SDB" "$RES" "$TMP/mm.tsv" --format-output "query,target,evalue,bits" >/dev/null
  awk 'BEGIN{FS=OFS="\t"}{k=$1; if(!(k in best) || $4+0>best[k]){best[k]=$4; tgt[k]=$2}} END{for(k in best) print k,tgt[k]}' \
    "$TMP/mm.tsv" > "$RECIP_FILE"
else
  : > "$RECIP_FILE"
fi

# ── 9) Emit final bridge table ───────────────────────────────────────────────
OUT_TSV="$OUTDIR/sarg_nt_bridge.tsv"
echo -e "prot_db\tprot_acc\tprot_len_aa\tprot_name\tnt_hash\tnt_len_bp\taa_identity\taa_aln_len\tqcov_aa\tscov_cds\tbitscore\tevalue\trank_status\tnear_top_hits\tdelta_bitscore_to_next\tframe\tstart_nt\tend_nt\t$RECIP_COL" > "$OUT_TSV"

# names & lengths
awk 'BEGIN{FS=">"; OFS="\t"} /^>/{sub(/^>/,""); hdr=$0; split(hdr,a," "); acc=a[1]; print acc"\t"hdr}' "$SARG_AA" > "$TMP/prot_headers.tsv"
seqkit fx2tab -n -l "$SARG_AA" | awk 'BEGIN{FS="\t"}{print $1"\t"$2}' > "$TMP/prot_lens.tsv"

awk -v RECIP="$RECIP" -v OFS="\t" '
  BEGIN{
    while((getline < "'"$TMP/prot_headers.tsv"'")>0){PNAME[$1]=$2}
    while((getline < "'"$TMP/prot_lens.tsv"'")>0){PLEN[$1]=$2}
    if(RECIP=="1"){ while((getline < "'"$RECIP_FILE"'")>0){RBQ=$1; RBS=$2; RB[RBQ]=RBS} }
  }
  # RANKED: qacc sseqid id_lenw cumAA qcov scov sumBits qlen slen qst qen sst sen status nb delta
  {
    q=$1; s=$2; pid=$3; cum=$4; qcov=$5; scov=$6; sb=$7; qlen=$8; slen=$9; qst=$10; qen=$11; sst=$12; sen=$13; status=$14; nb=$15; delta=$16
    frame=(sst<=sen)?"+":"-"
    recip_ok=(RECIP=="1" && RB[q]==s)?1:0
    print "SARG", q, (PLEN[q]+0), PNAME[q], s, slen, pid, cum, qcov, scov, sb, "NA", status, nb, delta, frame, sst, sen, recip_ok
  }
' "$RANKED" >> "$OUT_TSV"

# ── 10) Manifest ─────────────────────────────────────────────────────────────
{
  echo -e "file\tbytes\tsha256"
  for f in "$OUT_TSV" "$RAW_TSV" "$SORTED" "$AGG_TSV" "$FILTERED" "$RANKED"; do
    [[ -s "$f" ]] || continue
    sz=$(stat -c%s "$f"); sh=$(sha256sum "$f" | awk '{print $1}')
    echo -e "${f}\t${sz}\t${sh}"
  done
} > "$OUTDIR/manifest.tsv"

echo "[OK] Bridge table: $OUT_TSV"
echo "[OK] Rows: $(($(wc -l < "$OUT_TSV")-1))"
[[ "$RECIP" == "1" ]] && echo "[NOTE] Reciprocal flag via mmseqs (supportive only)."
