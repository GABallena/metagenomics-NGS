#!/usr/bin/env python3
"""Calculate plasmid-mapped read percentage.

Supports:
- Snakemake execution (snakemake.input / snakemake.output)
- Standalone CLI execution
"""

import argparse

def compute(plasmid_count_path, total_count_path, out_path):
    with open(plasmid_count_path) as f:
        plasmid_reads = int(f.read().strip())
    with open(total_count_path) as f:
        total_reads = int(f.read().strip())
    pct = (plasmid_reads / total_reads) * 100 if total_reads else 0.0
    with open(out_path, "w") as f:
        f.write(f"{pct:.2f}% of reads map to plasmids\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plasmid-read-count", required=True)
    parser.add_argument("--total-read-count", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    compute(args.plasmid_read_count, args.total_read_count, args.out)

if __name__ == "__main__":
    try:
        snakemake  # type: ignore
        compute(snakemake.input.plasmid_read_count, snakemake.input.total_read_count, snakemake.output.plasmid_percentage)
    except NameError:
        main()
