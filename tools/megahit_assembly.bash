#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Usage: bash megahit_assembly.bash <R1.fastq.gz> <R2.fastq.gz> <output_contigs.fasta>
if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> <output_contigs.fasta>" >&2
  exit 1
fi
# Input parameters
fasta1=$1
fasta2=$2
output_contigs=$3
output_dir=$(dirname $output_contigs)

# MEGAHIT settings
MEGAHIT_SETTINGS="--k-min 35 --k-max 141 --k-step 28"

# Ensure MEGAHIT output directory for this sample exists
mkdir -p "$output_dir"

# Remove any previous assembly directory for this sample
rm -rf "$output_dir/final_assembly"

# Run MEGAHIT
megahit -1 "$fasta1" -2 "$fasta2" -o "$output_dir/final_assembly" $MEGAHIT_SETTINGS --out-prefix final

# Check if MEGAHIT succeeded
if [ ! -f ""$output_dir/final_assembly/final.contigs.fa"" ]; then
  echo "Error: MEGAHIT failed or final.contigs.fa not found for inputs."
  exit 1
fi

# Filter contigs shorter than 200 bp
awk '/^>/ {if (seqlen >= 200) print seqname "\n" seq; seqname=$0; seqlen=0; seq=""; next} {seqlen += length($0); seq = seq $0} END {if (seqlen >= 200) print seqname "\n" seq}' \
""$output_dir/final_assembly/final.contigs.fa"" > "$output_contigs"

echo "$fasta1 and $fasta2 processing complete."
