#!/usr/bin/env python3
"""Extract dataset names from a BUSCO file_versions.tsv.

Public-safe helper utility.
"""

import argparse
import pandas as pd

def extract_all_datasets(input_file: str):
    df = pd.read_csv(input_file, sep="\t", header=None)
    df.columns = ["Dataset", "Date", "Hash", "Domain", "Type"]
    return df["Dataset"].tolist()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to file_versions.tsv")
    parser.add_argument("--output", default="all_dataset_list.txt")
    args = parser.parse_args()

    datasets = extract_all_datasets(args.input)
    with open(args.output, "w") as f:
        f.write("\n".join(datasets))

if __name__ == "__main__":
    main()
