#!/usr/bin/env python3
"""Compute contig quality metrics (N50/N90, GC, mean coverage from BAM).

Supports Snakemake and standalone CLI.
"""

import argparse
from Bio import SeqIO
import pysam

def calculate_nxx_lxx(lengths, frac):
    total = sum(lengths)
    target = total * frac
    run = 0
    lxx = 0
    nxx = 0
    for length in sorted(lengths, reverse=True):
        run += length
        lxx += 1
        if run >= target:
            nxx = length
            break
    return nxx, lxx

def gc_by_contig(records):
    out = {}
    for r in records:
        seq = str(r.seq).upper()
        out[r.id] = (seq.count("G") + seq.count("C")) / len(seq) * 100 if len(seq) else 0.0
    return out

def mean_coverage_by_contig(bam_path, contig_ids):
    cov = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for cid in contig_ids:
            total = 0
            covered = 0
            for pileupcolumn in bam.pileup(cid, truncate=True):
                total += pileupcolumn.n
                covered += 1
            cov[cid] = (total / covered) if covered else 0.0
    return cov

def run(contigs_fa, bam_file, quality_report, length_distribution):
    contigs = list(SeqIO.parse(contigs_fa, "fasta"))
    lengths = [len(c) for c in contigs]
    n50, l50 = calculate_nxx_lxx(lengths, 0.5)
    n90, l90 = calculate_nxx_lxx(lengths, 0.9)
    gc = gc_by_contig(contigs)
    cov = mean_coverage_by_contig(bam_file, [c.id for c in contigs])

    with open(quality_report, "w") as r:
        r.write(f"N50: {n50}\nL50: {l50}\nN90: {n90}\nL90: {l90}\n")
        r.write("Contig GC Content:\n")
        for c in contigs:
            r.write(f"{c.id}: {gc[c.id]:.2f}%\n")
        r.write("Contig Coverage:\n")
        for c in contigs:
            r.write(f"{c.id}: {cov.get(c.id, 0.0):.4f}x\n")

    with open(length_distribution, "w") as f:
        for length in sorted(lengths, reverse=True):
            f.write(f"{length}\n")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--contigs", required=True)
    p.add_argument("--bam", required=True)
    p.add_argument("--quality-report", required=True)
    p.add_argument("--length-distribution", required=True)
    args = p.parse_args()
    run(args.contigs, args.bam, args.quality_report, args.length_distribution)

if __name__ == "__main__":
    try:
        snakemake  # type: ignore
        run(snakemake.input.arg_contigs, snakemake.input.bam_file, snakemake.output.quality_report, snakemake.output.length_distribution)
    except NameError:
        main()
