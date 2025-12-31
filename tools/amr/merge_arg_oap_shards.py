#!/usr/bin/env python3
"""
Merge ARG-OAP Stage 2 per-shard outputs into per-sample files.

For each sample that has shard directories like samples/<SAMPLE>.part_XXX,
this script will find Stage 2 output tables (e.g., tpm.type.txt, rpkm.gene.txt,
normalized_cell.subtype.txt, etc.), sum values across shards per feature, and
write merged per-sample tables under samples/<SAMPLE>/merged/.

Assumptions
- Each shard table is a TSV with 3 columns: <feature> \t <R1-column> \t <R2-column>.
  The header names include "_R1_" and "_R2_" to disambiguate.
- Values are numeric; non-numeric values are treated as zero with a warning once per file.
- Summation is performed directly on the provided values. If you need to re-normalize
  (e.g., recompute TPM/RPKM from unnormalized counts), do that downstream.

Usage:
  python merge_arg_oap_shards.py --work-base arg_oap_work --samples SAMPLE-28_R13 SAMPLE-23_R8
  python merge_arg_oap_shards.py --work-base arg_oap_work   # auto-detect all samples with shards
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import OrderedDict
from glob import glob


METRICS = [
    "tpm",
    "rpkm",
    "ppm",
    "unnormalized_count",
    "unnormalized_copy",
    "normalized_cell",
    "normalized_16S",
]
LEVELS = ["gene", "subtype", "type"]


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def detect_samples_with_shards(samples_dir: str) -> list[str]:
    parts = [
        d
        for d in os.listdir(samples_dir)
        if os.path.isdir(os.path.join(samples_dir, d)) and ".part_" in d
    ]
    bases = sorted({re.sub(r"\.part_.*$", "", d) for d in parts})
    return bases


def find_shard_dirs(samples_dir: str, sample: str) -> list[str]:
    patt = os.path.join(samples_dir, f"{sample}.part_*")
    return sorted([d for d in glob(patt) if os.path.isdir(d)])


def list_present_tables(shards: list[str]) -> list[tuple[str, str, str]]:
    """Return list of (metric, level, filename) that exists in at least one shard.

    filename is like "tpm.type.txt".
    """
    present: set[tuple[str, str]] = set()
    for d in shards:
        for fn in os.listdir(d):
            for m in METRICS:
                for lvl in LEVELS:
                    if fn == f"{m}.{lvl}.txt":
                        present.add((m, lvl))
    out = [(m, lvl, f"{m}.{lvl}.txt") for (m, lvl) in sorted(present)]
    return out


def parse_numeric(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def merge_one_table(shards: list[str], filename: str, sample: str) -> tuple[list[str], list[list[str]]]:
    """Merge a single table over all shards.

    Returns (header, rows), where:
      header = [feature_name, f"{sample}_R1_trimmed", f"{sample}_R2_trimmed"]
      rows   = list of [feature_value, sum_R1, sum_R2] as strings
    """
    merged = OrderedDict()  # feature -> [sumR1, sumR2]
    feature_header = None
    warned_non_numeric = False

    for d in shards:
        path = os.path.join(d, filename)
        if not os.path.isfile(path):
            continue
        with open(path, "r", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration:
                continue
            if len(header) < 3:
                log(f"[WARN] Unexpected header with <3 columns in {path}: {header}")
                continue
            if feature_header is None:
                feature_header = header[0]
            # Identify R1 and R2 column positions
            # Expect exactly 2 numeric columns; use regex fallback
            r1_idx, r2_idx = 1, 2
            for i, h in enumerate(header[1:], start=1):
                if "_R1_" in h:
                    r1_idx = i
                if "_R2_" in h:
                    r2_idx = i

            for row in reader:
                if not row:
                    continue
                feat = row[0]
                v1 = parse_numeric(row[r1_idx]) if len(row) > r1_idx else float("nan")
                v2 = parse_numeric(row[r2_idx]) if len(row) > r2_idx else float("nan")
                if (v1 != v1) or (v2 != v2):  # NaN check
                    if not warned_non_numeric:
                        log(f"[WARN] Non-numeric encountered in {path}; treating as 0 for affected cells.")
                        warned_non_numeric = True
                    if v1 != v1:
                        v1 = 0.0
                    if v2 != v2:
                        v2 = 0.0
                if feat not in merged:
                    merged[feat] = [0.0, 0.0]
                merged[feat][0] += v1
                merged[feat][1] += v2

    if feature_header is None:
        # No shards had this table
        return [], []

    out_header = [feature_header, f"{sample}_R1_trimmed", f"{sample}_R2_trimmed"]
    out_rows: list[list[str]] = []
    for feat, (s1, s2) in merged.items():
        out_rows.append([feat, f"{s1}", f"{s2}"])

    return out_header, out_rows


def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def main() -> int:
    ap = argparse.ArgumentParser(description="Merge ARG-OAP per-shard outputs into per-sample files")
    ap.add_argument("--work-base", default="arg_oap_work", help="Path to ARG-OAP work dir (default: arg_oap_work)")
    ap.add_argument("--samples", nargs="*", help="Sample IDs to process (default: auto-detect by shards)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing merged files")
    args = ap.parse_args()

    samples_dir = os.path.join(args.work_base, "samples")
    if not os.path.isdir(samples_dir):
        log(f"[ERR ] Samples dir not found: {samples_dir}")
        return 2

    samples = args.samples or detect_samples_with_shards(samples_dir)
    if not samples:
        log("[WARN] No samples with shards detected. Nothing to do.")
        return 0

    for s in samples:
        log(f"== [{s}] merging per-shard tables ==")
        shards = find_shard_dirs(samples_dir, s)
        if not shards:
            log(f"[WARN] No shard dirs found for sample {s}; skipping.")
            continue
        present = list_present_tables(shards)
        if not present:
            log(f"[WARN] No recognized tables found in shards for {s}; skipping.")
            continue
        out_dir = os.path.join(samples_dir, s, "merged")
        ensure_dir(out_dir)

        for metric, level, fn in present:
            out_path = os.path.join(out_dir, fn)
            if os.path.exists(out_path) and not args.overwrite:
                log(f"  [SKIP] {fn} already exists")
                continue
            header, rows = merge_one_table(shards, fn, s)
            if not header:
                log(f"  [SKIP] {fn} not present in any shard")
                continue
            # Write output TSV
            with open(out_path, "w", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t")
                writer.writerow(header)
                writer.writerows(rows)
            log(f"  [OK  ] wrote {os.path.relpath(out_path)}")

        # Mark merge complete
        with open(os.path.join(out_dir, ".merged.done"), "w") as f:
            f.write("")

    log("[DONE] shard merging complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
