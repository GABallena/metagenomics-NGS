#!/usr/bin/env python3
"""Extract BUSCO lineage FASTA files named 'ancestral' and 'ancestral_variants'.

Public-safe helper utility.
"""

import os
import shutil
import argparse

def extract_fasta_files(source_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for root, _, files in os.walk(source_dir):
        for file in files:
            if file in ["ancestral", "ancestral_variants"]:
                lineage_name = os.path.basename(root)
                output_filename = f"{lineage_name}_{file}.fasta"
                shutil.copy(os.path.join(root, file), os.path.join(output_dir, output_filename))

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", default="busco_downloads/lineages")
    p.add_argument("--out", default="extracted_fastas")
    args = p.parse_args()
    extract_fasta_files(args.source, args.out)
    print("Done.")

if __name__ == "__main__":
    main()
