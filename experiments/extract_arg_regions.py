#!/usr/bin/env python3
"""Extract regions around ARG hits from BLAST tabular results.

Public-safe helper utility.
Assumes BLAST outfmt includes qseqid/sseqid and start/end columns; adjust indices as needed.
"""

from Bio import SeqIO
import pandas as pd
import argparse

def merge_regions(regions):
    regions.sort()
    merged = []
    for contig, start, end in regions:
        if not merged or merged[-1][0] != contig or merged[-1][2] < start:
            merged.append([contig, start, end])
        else:
            merged[-1][2] = max(merged[-1][2], end)
    return [(c, s, e) for c, s, e in merged]

def extract_arg_regions(contigs_file, blast_file, output_bed, output_fasta, flank=5000):
    """Extract merged regions ±flank around ARG hits."""
    blast_df = pd.read_csv(blast_file, sep='\t', header=None)

    regions = []
    for _, row in blast_df.iterrows():
        contig = row[0]
        # Default assumes start/end in columns 6 and 7 (0-based indexing).
        start = max(0, int(row[6]) - flank)
        end = int(row[7]) + flank
        if start > end:
            start, end = end, start
        regions.append((contig, start, end))

    regions = merge_regions(regions)

    with open(output_bed, 'w') as bed:
        for contig, start, end in regions:
            bed.write(f"{contig}\t{start}\t{end}\n")

    contigs = SeqIO.to_dict(SeqIO.parse(contigs_file, "fasta"))
    with open(output_fasta, 'w') as fasta:
        for contig, start, end in regions:
            if contig in contigs:
                seq = contigs[contig].seq[start:end]
                fasta.write(f">{contig}_{start}_{end}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--blast", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--flank", type=int, default=5000)
    args = parser.parse_args()

    extract_arg_regions(args.contigs, args.blast, args.output_bed, args.output_fasta, flank=args.flank)

if __name__ == "__main__":
    main()
