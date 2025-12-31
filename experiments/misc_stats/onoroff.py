#!/usr/bin/env python3
"""Classify genes as ON/OFF from a DESeq2-style results table.

Input: TSV with at least a log2 fold-change column (default: log2FoldChange)
Output: A small text summary + an annotated TSV

This is a simple, portfolio-safe utility script.
"""

from __future__ import annotations

import argparse
import pandas as pd
import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser(description="Mark genes as 'on'/'off' using a log2FC threshold.")
    ap.add_argument("--input", required=True, help="DESeq2 results TSV")
    ap.add_argument("--out-tsv", required=True, help="Output TSV with status column")
    ap.add_argument("--out-summary", required=True, help="Output summary text file")
    ap.add_argument("--log2fc-col", default="log2FoldChange", help="Column name for log2FC (default: log2FoldChange)")
    ap.add_argument("--threshold", type=float, default=-1.0, help="OFF threshold for log2FC (default: -1.0)")
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    if args.log2fc_col not in df.columns:
        raise ValueError(f"Column '{args.log2fc_col}' not found. Available: {list(df.columns)}")

    df["status"] = np.where(df[args.log2fc_col] <= args.threshold, "off", "on")

    total = len(df)
    off = int((df["status"] == "off").sum())
    pct_off = (off / total * 100.0) if total else 0.0

    df.to_csv(args.out_tsv, sep="\t", index=False)

    with open(args.out_summary, "w", encoding="utf-8") as f:
        f.write(f"threshold\t{args.threshold}\n")
        f.write(f"total_genes\t{total}\n")
        f.write(f"off_genes\t{off}\n")
        f.write(f"percent_off\t{pct_off:.2f}\n")

    print(f"Wrote: {args.out_tsv}")
    print(f"Wrote: {args.out_summary}")


if __name__ == "__main__":
    main()
