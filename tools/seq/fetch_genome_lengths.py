#!/usr/bin/env python3
"""
Fetch genome lengths (total sequence length) per NCBI taxid using NCBI Datasets CLI.

Outputs a TSV: taxid, genome_length, ploidy

Usage examples:
  - From merged TSV produced by merge_metaphlan_bracken.py:
      python3 fetch_genome_lengths.py --merged diversity_results/metaphlan_bracken_merged.tsv --out genome-lengths.tsv
  - Directly from Bracken outputs:
      python3 fetch_genome_lengths.py --bracken-dir bracken_output --rank species --out genome-lengths.tsv

Requirements:
  - NCBI Datasets CLI available in PATH (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)

Selection policy for assemblies (per taxid):
  1) Prefer RefSeq assemblies, category Reference or Representative
  2) Prefer higher assembly_level: Complete Genome > Chromosome > Scaffold > Contig
  3) Prefer latest (by release_date or last_updated)
  4) Fallback to any with assembly_stats.total_sequence_length

Notes:
  - Ploidy is set to 1.0 unless you provide your own ploidy mapping.
  - For taxa with multiple species-level taxids, supply species taxids (rank=species recommended).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from typing import Dict, Iterable, List, Optional, Tuple

ASSEMBLY_LEVEL_ORDER = {
    "Complete Genome": 4,
    "Chromosome": 3,
    "Scaffold": 2,
    "Contig": 1,
}


def have_datasets_cli() -> bool:
    try:
        subprocess.run(["datasets", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        return True
    except Exception:
        return False


def pick_best_assembly(entries: List[dict]) -> Optional[dict]:
    if not entries:
        return None
    # Score and pick
    best = None
    # Tuple: (has_len:int, cat_ref_bonus:int, level_score:int, updated:str)
    best_score = (-1, -1, -1, "")
    for e in entries:
        asm = e.get("assembly", {})
        stats = e.get("assembly_stats", {})
        org = e.get("organism", {})
        source_db = asm.get("assembly_source") or asm.get("assembly_db") or ""
        cat = asm.get("assembly_category") or ""
        level = asm.get("assembly_level") or ""
        level_score = ASSEMBLY_LEVEL_ORDER.get(level, 0)
        # Prefer RefSeq (or genbank if no refseq)
        refseq_bonus = 1 if str(source_db).lower() == "refseq" else 0
        # Prefer reference/representative
        cat_bonus = 2 if cat in ("reference", "representative") else 0
        # Prefer those that actually have total length
        length_raw = stats.get("total_sequence_length")
        try:
            length = float(length_raw) if length_raw is not None else 0.0
        except Exception:
            length = 0.0
        has_len = 1 if length > 0 else 0
        # Tie-breaker by last-updated string
        updated = asm.get("last_updated") or asm.get("submission_date") or ""
        score = (int(has_len), int(cat_bonus + refseq_bonus), int(level_score), str(updated or ""))
        if score > best_score:
            best_score = score
            best = e
    return best


def get_genome_length_for_taxid(taxid: str) -> Optional[int]:
    """Query NCBI Datasets CLI for a taxid and return total_sequence_length of the best assembly."""
    try:
        # Use JSON lines to stream results
        proc = subprocess.run(
            ["datasets", "summary", "genome", "taxon", str(taxid), "--as-json-lines"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            sys.stderr.write(f"WARN: datasets summary failed for taxid {taxid}: {proc.stderr.strip()}\n")
            return None
        entries: List[dict] = []
        for line in proc.stdout.splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
                # Depending on version, data may be under 'assemblies' or top-level
                if isinstance(obj, dict) and "assemblies" in obj:
                    for a in obj.get("assemblies", []) or []:
                        entries.append(a)
                else:
                    entries.append(obj)
            except Exception:
                continue
        best = pick_best_assembly(entries)
        if not best:
            return None
        stats = best.get("assembly_stats", {})
        total = stats.get("total_sequence_length")
        if isinstance(total, (int, float)):
            return int(total)
        try:
            return int(float(total)) if total is not None else None
        except Exception:
            return None
    except Exception as e:
        sys.stderr.write(f"WARN: failed datasets query for taxid {taxid}: {e}\n")
        return None


def collect_taxids_from_merged(path: str, rank: str) -> List[str]:
    taxids: List[str] = []
    rank = rank.lower()
    with open(path, "r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if (row.get("rank", "").lower() == rank) and row.get("taxid"):
                taxids.append(str(row["taxid"]).strip())
    return sorted({t for t in taxids if t and t != "nan"})


def collect_taxids_from_bracken(bracken_dir: str, rank: str) -> List[str]:
    taxids: List[str] = []
    rank = rank.lower()
    suffix = f".{rank}.bracken"
    breport = f".{rank}.breport"
    for root, _dirs, files in os.walk(bracken_dir):
        for name in files:
            low = name.lower()
            if low.endswith(suffix) or low.endswith(breport):
                p = os.path.join(root, name)
                try:
                    with open(p, "r") as fh:
                        for line in fh:
                            if not line or line.startswith("#"):
                                continue
                            cols = line.rstrip("\n").split("\t")
                            if len(cols) >= 7:
                                # .bracken format: name taxid lvl ...
                                taxids.append(cols[1].strip())
                            elif len(cols) == 6:
                                # .breport: % clade_reads taxon_reads rank taxid name
                                taxids.append(cols[4].strip())
                except Exception:
                    pass
    return sorted({t for t in taxids if t and t != "nan"})


def main():
    ap = argparse.ArgumentParser(description="Fetch genome lengths from NCBI Datasets for taxids.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--merged", help="Path to merged TSV from merge_metaphlan_bracken.py")
    src.add_argument("--bracken-dir", help="Path to directory with Bracken outputs")
    ap.add_argument("--rank", default="species", choices=["species","genus","family","order","class","phylum"], help="Rank to target when collecting taxids (default species)")
    ap.add_argument("--out", default="genome-lengths.tsv", help="Output TSV (taxid, genome_length, ploidy)")
    ap.add_argument("--limit", type=int, default=None, help="Optional limit of taxids to query (for testing)")
    args = ap.parse_args()

    if not have_datasets_cli():
        sys.stderr.write("ERROR: 'datasets' CLI not found in PATH. Install NCBI Datasets CLI and retry.\n")
        sys.exit(2)

    if args.merged:
        if not os.path.exists(args.merged):
            sys.stderr.write(f"ERROR: merged file not found: {args.merged}\n")
            sys.exit(2)
        taxids = collect_taxids_from_merged(args.merged, args.rank)
    else:
        if not os.path.isdir(args.bracken_dir):
            sys.stderr.write(f"ERROR: bracken-dir not found: {args.bracken_dir}\n")
            sys.exit(2)
        taxids = collect_taxids_from_bracken(args.bracken_dir, args.rank)

    if args.limit:
        taxids = taxids[: args.limit]

    if not taxids:
        sys.stderr.write("WARN: no taxids found to query.\n")

    out_path = args.out
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["taxid", "genome_length", "ploidy"])
        for tx in taxids:
            length = get_genome_length_for_taxid(tx)
            if length is None:
                sys.stderr.write(f"WARN: no genome length for taxid {tx}\n")
                continue
            w.writerow([tx, int(length), 1.0])
            fh.flush()
            sys.stdout.write(f"OK: {tx}\t{length}\n")

    print(f"Wrote: {out_path} (rows={sum(1 for _ in open(out_path)) - 1})")


if __name__ == "__main__":
    main()
