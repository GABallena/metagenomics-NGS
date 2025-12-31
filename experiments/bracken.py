#!/usr/bin/env python3
"""Run Bracken on a directory of Kraken2 .k2report files.

Public-safe wrapper: paths and parameters are provided via CLI args.
"""

import os
import subprocess
import argparse

def run_bracken(kraken_report_dir, bracken_out_dir, kraken_db, read_length, level, threads):
    os.makedirs(bracken_out_dir, exist_ok=True)

    reports = [f for f in os.listdir(kraken_report_dir) if f.endswith(".k2report")]
    if not reports:
        print("No .k2report files found. Run Kraken2 first.")
        return 0

    for report in reports:
        sample = report.replace(".k2report", "")
        bracken_input = os.path.join(kraken_report_dir, report)
        bracken_output = os.path.join(bracken_out_dir, f"{sample}.bracken")
        breport_output = os.path.join(bracken_out_dir, f"{sample}.breport")

        cmd = [
            "bracken",
            "-d", kraken_db,
            "-i", bracken_input,
            "-o", bracken_output,
            "-w", breport_output,
            "-r", str(read_length),
            "-l", level,
            "-t", str(threads)
        ]

        print(f"Running Bracken for {sample}...")
        subprocess.run(cmd, check=True)

    print("Bracken completed.")
    return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kraken-report-dir", default="Kraken/kraken_output")
    parser.add_argument("--bracken-out-dir", default="Kraken/bracken_output")
    parser.add_argument("--kraken-db", default="Kraken/k2_db")
    parser.add_argument("--read-length", type=int, default=100)
    parser.add_argument("--level", default="S")
    parser.add_argument("--threads", type=int, default=10)
    args = parser.parse_args()
    raise SystemExit(run_bracken(args.kraken_report_dir, args.bracken_out_dir, args.kraken_db, args.read_length, args.level, args.threads))

if __name__ == "__main__":
    main()
