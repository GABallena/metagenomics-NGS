#!/usr/bin/env bash
#
# Portfolio-safe copy: paths/identifiers generalized; databases not included.
#
set -euo pipefail

IN=work/argdb_harmonized/nonredundant.fa
OUT=work/argdb_harmonized/nonredundant.clean.fa

awk '
  BEGIN{FS=OFS=""}
  /^>/ {print; blank_ok=1; next}
  NF==0 && blank_ok {next}        # drop the blank line right after a header
  { gsub(/[ \t\r]/,""); if(NF){print; blank_ok=0} }
' "$IN" > "$OUT"

CLSTR=work/argdb_harmonized/cdhit_100.clstr
OUT=work/argdb_harmonized/cluster_map.fixed.tsv

# write the awk program to a file (prevents quote/escape issues)
cat > /tmp/rebuild_cluster_map.awk <<'AWK'
BEGIN{
  OFS="\t";
  print "orig_hash","rep_hash";
  rep="";
  n=0;
}
# Start of a new cluster
/^>Cluster/{
  if(n>0){
    if(rep=="") rep=members[1];
    for(i=1;i<=n;i++){ print members[i], rep }
  }
  rep="";
  n=0;
  next;
}
# Member lines, e.g.:
# "0    891nt, >fa1903... at 1:891... *"
{
  if (match($0, />[0-9a-f]{64}/)) {
    hash=substr($0, RSTART+1, RLENGTH-1);
    members[++n]=hash;
    if(index($0,"*")>0){ rep=hash }
  }
}
END{
  if(n>0){
    if(rep=="") rep=members[1];
    for(i=1;i<=n;i++){ print members[i], rep }
  }
}
AWK

# run it
awk -f /tmp/rebuild_cluster_map.awk "$CLSTR" > "$OUT"
FA=work/argdb_harmonized/nonredundant.clean.fa
[[ -s "$FA" ]] || FA=work/argdb_harmonized/nonredundant.fa

# unique FASTA ids
grep '^>' "$FA" | sed 's/^>//' | sort -u > /tmp/fa.ids
# unique rep ids in fixed map
cut -f2 "$OUT" | awk 'NR>1' | sort -u > /tmp/rep.ids

echo "rep_hash not present in FASTA (should be 0):"
comm -23 /tmp/rep.ids /tmp/fa.ids | wc -l
comm -23 /tmp/rep.ids /tmp/fa.ids | head -n 10

echo "FASTA ids not listed as rep_hash (info only):"
comm -23 /tmp/fa.ids /tmp/rep.ids | wc -l | sed 's/$/ (ok if >0)/'
mv work/argdb_harmonized/cluster_map.fixed.tsv work/argdb_harmonized/cluster_map.tsv

