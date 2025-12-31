#!/usr/bin/env python3
"""
Scaffolding utility for extending contigs using long reads.

Assumptions:
- Runs in a controlled environment where required tools are available
- Input/output directories are relative to the project root
- Not intended to be executed with elevated privileges on shared HPC systems
"""

import os
import subprocess
from pathlib import Path

LONG_READS_DIR = Path("long_reads")
FILTERED_CONTIGS_DIR = Path("filtered_contigs")
SCAFFOLDS_DIR = Path("scaffolds")
LOGS_DIR = Path("logs")

for d in [SCAFFOLDS_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

def run_command(cmd, log_file):
    with open(log_file, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=log, check=True)

for contig_file in FILTERED_CONTIGS_DIR.glob("*.fasta"):
    sample = contig_file.stem
    long_reads = LONG_READS_DIR / f"{sample}.fastq.gz"
    output_fasta = SCAFFOLDS_DIR / f"{sample}_scaffolded.fasta"
    log_file = LOGS_DIR / f"{sample}.log"

    cmd = [
        "contig-extender/dist/extender_wrapper",
        "-c", str(contig_file),
        "-r", str(long_reads),
        "-o", str(output_fasta)
    ]

    run_command(cmd, log_file)

print("Scaffolding completed.")
