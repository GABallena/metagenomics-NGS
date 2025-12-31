#!/usr/bin/env python3
"""Move final filtered contigs from per-sample assembly outputs into a single directory.

Public-safe utility: uses relative defaults and accepts CLI arguments.
"""

import os
import shutil
import glob
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-dir", default="spades_output", help="Directory containing per-sample assembly folders")
    parser.add_argument("--target-dir", default="filtered_contigs", help="Destination directory for contig FASTA files")
    parser.add_argument("--pattern", default="*_assembly/*_final_filtered_contigs.fa", help="Glob pattern relative to assembly-dir")
    args = parser.parse_args()

    os.makedirs(args.target_dir, exist_ok=True)

    contig_files = glob.glob(os.path.join(args.assembly_dir, args.pattern))
    if not contig_files:
        print("No contig files matched. Check --assembly-dir and --pattern.")
        return 0

    for file_path in contig_files:
        file_name = os.path.basename(file_path)
        dest_path = os.path.join(args.target_dir, file_name)
        shutil.move(file_path, dest_path)
        print(f"Moved: {file_path} -> {dest_path}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
