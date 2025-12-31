#!/bin/bash
#
# Portfolio-safe copy: paths/identifiers generalized; databases not included.
#
set -euo pipefail

MAX_RECORDS=1000
OUTDIR="ncbi_mge_fastas"
mkdir -p "$OUTDIR"

keywords=(
"IS26" "ISCR1" "Tn3" "Tn21" "Tn7" "Tn4401" "Tn125"
"IncF" "IncA/C2" "IncX3" "ISAba1" "intI1" "Tn6029"
"IS6100" "IS1R" "ISKpn19" "ISPa1328" "Tn2" "IncHI2"
"IS5075" "ISCR2" "Tn1721" "Tn3-family" "IS4321" "intI2"
"IncN" "IncI1" "ISCR4" "ISCR10" "ISKpn7" "IS2" "IS3"
"IS911" "IS4" "IncQ1" "intI3" "Tn10" "IS407" "Tn554"
"IncL/M" "IS5" "IS10" "IS1414-like" "ISaz1-like" "ISBha1"
"ISGpn1" "Desulfotalea-like IS" "IS605-like cryptic"
"RSF1010-type replication genes" "MITEs" "Foldback elements"
"Truncated IS remnants" "Phage tail/capsid genes only"
"Partial tnpA-like ORFs" "RecA / RuvC-like RNase genes"
"RelBE/Phd-Doc TA systems" "Type II RM system genes"
"Alu-like or IS-like repeats" "Chromosomal integrase-like ORFs"
)

echo "[*] Starting NCBI sequence fetch for expanded MGE keyword list..."

for keyword in "${keywords[@]}"; do
    clean_name=$(echo "$keyword" | tr -s ' /-' '__' | tr -cd '[:alnum:]_')
    echo "→ Searching for: '$keyword' → Output: $OUTDIR/${clean_name}.fasta"

    esearch -db nuccore -query "$keyword[All Fields]" |
        efetch -format fasta |
        head -n $((MAX_RECORDS * 2)) > "$OUTDIR/${clean_name}.fasta"

    sleep 1
    if [[ ! -s "$OUTDIR/${clean_name}.fasta" ]]; then
        echo "⚠️  Warning: No sequences found for '$keyword'."
    else
        echo "[✓] Saved $(grep -c '^>' "$OUTDIR/${clean_name}.fasta") sequences for '$keyword'."
    fi
done

echo "[✓] All sequences fetched into: $OUTDIR/"
