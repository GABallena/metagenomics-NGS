#!/usr/bin/env python3
"""
Merge MetaPhlAn and Bracken outputs into a single table.

Inputs (defaults assume this repo layout):
- MetaPhlAn profiles: metaphlan_out/*_profile.txt (also supports .biom/.json BIOM)
- Bracken outputs at chosen rank: bracken_output/*.{rank}.bracken (or .breport)

Outputs:
- metaphlan_bracken_merged.tsv (long/tidy): sample, taxon, rank, metaphlan_abundance, bracken_abundance
- metaphlan_bracken_merged_wide.tsv (wide):  sample, taxon, rank, metaphlan, bracken

Notes:
- Abundances are percentages for both tools (MetaPhlAn relative_abundance; Bracken fraction_total_reads*100).
- Taxon names are normalized to plain names (e.g., "Escherichia coli").
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np


# Map rank keyword to MetaPhlAn/Bracken conventions
RANK_PREFIX = {
    "species": "s__",
    "genus": "g__",
    "family": "f__",
    "order": "o__",
    "class": "c__",
    "phylum": "p__",
    "kingdom": "k__",
    "domain": "d__",
}

BRACKEN_RANK_SUFFIX = {
    "species": "species",
    "genus": "genus",
    "family": "family",
    "order": "order",
    "class": "class",
    "phylum": "phylum",
}

BRACKEN_RANK_LETTER = {
    "species": "S",
    "genus": "G",
    "family": "F",
    "order": "O",
    "class": "C",
    "phylum": "P",
    "kingdom": "K",
    "domain": "D",
}

# Canonical rank order for nicer sorting (polish)
RANK_ORDER = ["phylum", "class", "order", "family", "genus", "species"]

# Parsing constants (document magic numbers)
BRACKEN_MIN_COLUMNS = 7
METAPHLAN_ABUNDANCE_FALLBACK_INDEX = 2
def detect_metaphlan_header(header_line: str) -> Optional[List[str]]:
    """Detect and validate MetaPhlAn header format; return columns if valid."""
    parts = header_line.lstrip("# ").split("\t")
    required_candidates = {"clade_name", "relative_abundance"}
    if any(col in parts for col in required_candidates):
        return parts
    return None


def canonicalize_taxon_name(name: str) -> str:
    """Light canonization to align names across tools.

    - Replace underscores with spaces
    - Remove bracketed/parenthetical annotations (possible strain info)
    - Drop 'Candidatus'/'Ca.' prefixes
    - Strip quotes and collapse whitespace
    """
    if not name:
        return name
    s = name.replace("_", " ")
    s = re.sub(r"\s*[\(\[].*?[\)\]]\s*", " ", s)
    s = re.sub(r"^(?:Candidatus\s+|Ca\.\s+)", "", s)
    s = s.replace('"', "").strip()
    s = re.sub(r"\s+", " ", s)
    return s

def infer_sample_id(path: str) -> str:
    """Get sample ID from file name like SAMPLE-01_R1_profile.txt or SAMPLE-01_R1.species.bracken."""
    base = os.path.basename(path)
    # Strip known suffixes
    base = re.sub(r"_profile\.txt$", "", base)
    base = re.sub(r"\.(species|genus|family|order|class|phylum)\.bracken$", "", base)
    base = re.sub(r"\.(species|genus|family|order|class|phylum)\.breport$", "", base)
    base = re.sub(r"\.biom$", "", base, flags=re.IGNORECASE)
    base = re.sub(r"\.json$", "", base, flags=re.IGNORECASE)
    return base

def normalize_taxon_from_metaphlan(clade_name: str, rank: str) -> Optional[str]:
    """Extract plain taxon name from MetaPhlAn clade_name at the requested rank.

    clade_name examples: "k__Bacteria|p__Firmicutes|c__Bacilli|...|s__Escherichia_coli"
    Returns: "Escherichia coli" or None if the clade doesn't contain the rank.
    """
    prefix = RANK_PREFIX.get(rank)
    if not prefix:
        return None
    parts = clade_name.split("|")
    # Find last part matching prefix, in case multiple present
    for part in reversed(parts):
        if part.startswith(prefix):
            name = part[len(prefix) :]
            # Ignore unclassified/unknown
            if not name or name.lower() in {"unclassified", "unknown", "nan"}:
                return None
            return name.replace("_", " ")
    return None

def _is_json_biom(path: str) -> bool:
    """Quickly detect if a file is BIOM JSON (starts with { or [ after whitespace)."""
    try:
        with open(path, "rb") as f:
            head = f.read(8192)  # read more to tolerate leading whitespace/newlines
            c = head.lstrip()
            return c[:1] in (b"{", b"[")
    except Exception:
        return False

def _parse_metaphlan_biom_json(profile_path: str, rank: str) -> pd.DataFrame:
    """Parse a MetaPhlAn BIOM JSON profile into a DataFrame.

    Returns columns: [sample, taxon, rank, metaphlan_abundance, taxid]
    Note: JSON profile rows typically lack NCBI taxids; taxid will be None.
    """
    import json

    with open(profile_path, "rb") as f:
        j = json.load(f)

    # Sample id from BIOM columns, fallback to filename
    sample = None
    try:
        cols = j.get("columns") or []
        if cols and isinstance(cols, list) and isinstance(cols[0], dict):
            sample = cols[0].get("id") or cols[0].get("sample_id")
    except Exception:
        sample = None
    if not sample:
        sample = infer_sample_id(profile_path)

    rows = j.get("rows") or []
    data = j.get("data") or []
    matrix_type = j.get("matrix_type", "sparse")

    records: List[Tuple[str, str, str, float]] = []

    if matrix_type == "sparse":
        # data as list of [row_idx, col_idx, value]
        for triplet in data:
            try:
                ri, ci, val = triplet
            except Exception:
                continue
            # Only one column expected; ignore ci
            try:
                clade = rows[ri].get("id")
            except Exception:
                continue
            taxon = normalize_taxon_from_metaphlan(clade, rank)
            if not taxon:
                continue
            try:
                abundance = float(val)
            except Exception:
                continue
            records.append((sample, canonicalize_taxon_name(taxon), rank, abundance, None))
    else:
        # Dense matrix: expect data as 2D array of shape [n_rows, n_cols]
        try:
            for ri, rowvals in enumerate(data):
                if not rowvals:
                    continue
                try:
                    clade = rows[ri].get("id")
                except Exception:
                    continue
                taxon = normalize_taxon_from_metaphlan(clade, rank)
                if not taxon:
                    continue
                try:
                    abundance = float(rowvals[0])
                except Exception:
                    continue
                records.append((sample, canonicalize_taxon_name(taxon), rank, abundance, None))
        except Exception:
            pass

    if not records:
        return pd.DataFrame(columns=["sample", "taxon", "rank", "metaphlan_abundance", "taxid"])
    return pd.DataFrame(records, columns=["sample", "taxon", "rank", "metaphlan_abundance", "taxid"])

def parse_metaphlan_profile(profile_path: str, rank: str) -> pd.DataFrame:
    """Parse a single MetaPhlAn profile.txt into a DataFrame with columns:
    [sample, taxon, rank, metaphlan_abundance]
    """
    # If this is a BIOM JSON file (common in MetaPhlAn outputs saved as .biom or misnamed), parse accordingly
    if _is_json_biom(profile_path):
        return _parse_metaphlan_biom_json(profile_path, rank)

    sample = infer_sample_id(profile_path)
    records: List[Tuple[str, str, str, float]] = []
    with open(profile_path, "r") as fh:
        header_cols: Optional[List[str]] = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#") and header_cols is None:
                cols = detect_metaphlan_header(line)
                if cols is not None:
                    header_cols = cols
                continue
            # Split row manually to avoid csv field limits
            cols = line.split("\t")
            if header_cols is None:
                # Fallback minimal header
                header_cols = ["clade_name", "NCBI_tax_id", "relative_abundance"]
            # Map first len(header) columns
            entry: Dict[str, str] = {}
            for i in range(min(len(header_cols), len(cols))):
                entry[header_cols[i]] = cols[i]
            clade = entry.get("clade_name") or entry.get("clade") or (cols[0] if cols else None)
            taxid = entry.get("NCBI_tax_id") or (cols[1] if len(cols) > 1 else None)
            ra_str = entry.get("relative_abundance") or entry.get("relative_abundance_percentage")
            if ra_str is None and len(cols) > METAPHLAN_ABUNDANCE_FALLBACK_INDEX:
                ra_str = cols[METAPHLAN_ABUNDANCE_FALLBACK_INDEX]
            if clade is None or ra_str is None:
                continue
            taxon = normalize_taxon_from_metaphlan(clade, rank)
            if not taxon:
                continue
            try:
                ra = float(ra_str)
            except ValueError:
                continue
            records.append((sample, taxon, rank, ra, taxid))
    if not records:
        return pd.DataFrame(columns=["sample", "taxon", "rank", "metaphlan_abundance", "taxid"])
    df = pd.DataFrame(records, columns=["sample", "taxon", "rank", "metaphlan_abundance", "taxid"])
    # Canonicalize taxon naming
    df["taxon"] = df["taxon"].map(canonicalize_taxon_name)
    return df

def collect_bacterial_taxa_from_metaphlan(profile_path: str, ranks: List[str]) -> Dict[str, set]:
    """Scan a MetaPhlAn profile and collect taxa names at given ranks that belong to Bacteria.

    Returns dict rank -> set of normalized taxon names.
    """
    allowed: Dict[str, set] = {rk: set() for rk in ranks}
    # JSON BIOM support
    if _is_json_biom(profile_path):
        import json
        try:
            with open(profile_path, "rb") as f:
                j = json.load(f)
            rows = j.get("rows") or []
            for row in rows:
                clade = row.get("id") or ""
                # Accept either kingdom or domain prefix for Bacteria
                if ("k__Bacteria" not in clade) and ("d__Bacteria" not in clade):
                    continue
                for rk in ranks:
                    taxon = normalize_taxon_from_metaphlan(clade, rk)
                    if taxon:
                        allowed[rk].add(canonicalize_taxon_name(taxon))
        except Exception as e:
            print(f"WARN: failed to scan JSON BIOM {profile_path}: {e}")
        return allowed

    # Fallback: text profile
    with open(profile_path, "r") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if not cols:
                continue
            clade = cols[0]
            # Accept either kingdom or domain prefix for Bacteria
            if ("k__Bacteria" not in clade) and ("d__Bacteria" not in clade):
                continue
            for rk in ranks:
                taxon = normalize_taxon_from_metaphlan(clade, rk)
                if taxon:
                    allowed[rk].add(canonicalize_taxon_name(taxon))
    return allowed

def collect_bacterial_taxa(meta_files: List[str], ranks: List[str]) -> Dict[str, set]:
    """Aggregate bacterial taxa across multiple MetaPhlAn profiles for given ranks."""
    agg: Dict[str, set] = {rk: set() for rk in ranks}
    for path in meta_files:
        try:
            res = collect_bacterial_taxa_from_metaphlan(path, ranks)
            for rk in ranks:
                agg[rk].update(res.get(rk, set()))
        except Exception as e:
            print(f"WARN: failed to scan {path} for bacterial taxa: {e}")
    return agg

def parse_bracken_file(bracken_path: str, rank: str) -> pd.DataFrame:
    """Parse a single Bracken file into DataFrame with columns:
    [sample, taxon, rank, bracken_abundance]

    Expects columns: name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads
    """
    sample = infer_sample_id(bracken_path)
    # Determine rank from file suffix if possible
    m = re.search(r"\.(species|genus|family|order|class|phylum)\.(bracken|breport)$", bracken_path, re.IGNORECASE)
    file_rank = m.group(1) if m else None
    expected_letter = BRACKEN_RANK_LETTER.get(rank)

    records: List[Tuple[str, str, str, float]] = []
    with open(bracken_path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            # Defensive: some bracken outputs include a header; handle both
            if (row[0] == "name" and len(row) >= 7) or (row[0].lower() in {"%", "percent"} and len(row) >= 6):
                continue
            name = None
            taxid = None
            lvl = None
            abundance = None
            taxon = None
            if len(row) >= 7:
                # .bracken format
                name, taxid, lvl = row[0], row[1], row[2]
                frac_str = row[-1]
                try:
                    frac = float(frac_str)
                except ValueError:
                    continue
                abundance = frac * 100.0
                taxon = canonicalize_taxon_name((name or "").strip())
            elif len(row) == 6:
                # .breport format: %, clade_reads, taxon_reads, rank, taxid, name
                pct_str, _clade_reads, _taxon_reads, lvl, taxid, name = row
                try:
                    pct = float(pct_str)
                except ValueError:
                    continue
                abundance = pct  # already percent
                taxon = canonicalize_taxon_name((name or "").strip())
            else:
                continue
            # Filter by rank (case-insensitive in file)
            if expected_letter and (lvl or "").upper() != expected_letter:
                continue
            if not name or not taxon:
                continue
            # Basic cleanup: avoid generic labels
            if taxon.lower() in {"unclassified", "root", "no rank"}:
                continue
            records.append((sample, taxon, rank, abundance, taxid))
    if not records:
        return pd.DataFrame(columns=["sample", "taxon", "rank", "bracken_abundance", "taxid"])
    df = pd.DataFrame(records, columns=["sample", "taxon", "rank", "bracken_abundance", "taxid"])
    # Warn if file rank doesn't match requested rank
    if file_rank and file_rank != rank:
        print(f"WARN: {os.path.basename(bracken_path)} is '{file_rank}', but --rank is '{rank}'.")
    return df

def merge_metaphlan_bracken(metaphlan_dir: str, bracken_dir: str, rank: str, sample_glob: Optional[str] = None, join_on: str = "name") -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Produce merged long and wide DataFrames for the chosen rank.

    sample_glob can restrict which samples to include (pattern applied to filenames).
    """
    # Collect MetaPhlAn profiles: support text profiles and BIOM JSON saved as .biom/.json
    meta_pattern_txt = os.path.join(metaphlan_dir, "*_profile.txt")
    meta_pattern_biom = os.path.join(metaphlan_dir, "*.biom")
    meta_pattern_json = os.path.join(metaphlan_dir, "*.json")
    meta_files = sorted(set(glob.glob(meta_pattern_txt) + glob.glob(meta_pattern_biom) + glob.glob(meta_pattern_json)))
    if sample_glob:
        meta_files = [p for p in meta_files if re.search(sample_glob, os.path.basename(p))]

    # Collect Bracken files for rank (.bracken and .breport) and choose one per sample (prefer .bracken)
    br_suffix = BRACKEN_RANK_SUFFIX.get(rank, rank)
    br_pattern1 = os.path.join(bracken_dir, f"*.{br_suffix}.bracken")
    br_pattern2 = os.path.join(bracken_dir, f"*.{br_suffix}.breport")
    all_br = sorted(set(glob.glob(br_pattern1) + glob.glob(br_pattern2)))
    if sample_glob:
        all_br = [p for p in all_br if re.search(sample_glob, os.path.basename(p))]
    chosen: Dict[str, str] = {}
    seen_types: Dict[str, set] = {}
    for p in all_br:
        s = infer_sample_id(p)
        seen_types.setdefault(s, set()).add(os.path.splitext(p)[1].lower())
        if s not in chosen:
            chosen[s] = p
        else:
            # Prefer .bracken over .breport
            if p.lower().endswith('.bracken') and chosen[s].lower().endswith('.breport'):
                chosen[s] = p
    # Warn if both existed for any sample and .breport was ignored
    for s, types in seen_types.items():
        if {'.bracken', '.breport'}.issubset(types):
            print(f"WARN: Both .bracken and .breport found for sample {s} at rank {rank}; using .bracken")
    br_files = sorted(chosen.values())

    # Parse all
    meta_dfs: List[pd.DataFrame] = []
    for path in meta_files:
        try:
            meta_dfs.append(parse_metaphlan_profile(path, rank))
        except Exception as e:
            print(f"WARN: failed to parse MetaPhlAn file {path}: {e}")

    br_dfs: List[pd.DataFrame] = []
    for path in br_files:
        try:
            br_dfs.append(parse_bracken_file(path, rank))
        except Exception as e:
            print(f"WARN: failed to parse Bracken file {path}: {e}")

    # Filter out empty frames to avoid FutureWarning in pandas concat
    meta_dfs = [df for df in meta_dfs if df is not None and not df.empty]
    br_dfs = [df for df in br_dfs if df is not None and not df.empty]
    meta_df = pd.concat(meta_dfs, ignore_index=True) if meta_dfs else pd.DataFrame(columns=["sample","taxon","rank","metaphlan_abundance","taxid"]) 
    br_df = pd.concat(br_dfs, ignore_index=True) if br_dfs else pd.DataFrame(columns=["sample","taxon","rank","bracken_abundance","taxid"]) 

    # Taxon names already canonicalized in parsers

    # Merge
    if join_on == "taxid":
        # Warn if MetaPhlAn has no taxids (e.g., BIOM profiles)
        # This can lead to minimal or no matches when joining on taxid.
        try:
            if meta_df.shape[0] > 0 and meta_df["taxid"].isna().all():
                print("WARN: MetaPhlAn profiles contain no taxids; joining on taxid may yield few or no matches.")
        except Exception:
            pass
        # Merge by sample+rank+taxid; coalesce taxon names after merge
        merged = pd.merge(
            meta_df,
            br_df,
            on=["sample", "rank", "taxid"],
            how="outer",
            suffixes=("_m", "_b"),
        )
        # Choose a representative taxon name
        merged["taxon"] = merged.get("taxon_m").fillna(merged.get("taxon_b"))
        # Align abundance columns to canonical names
        if "metaphlan_abundance_m" in merged.columns:
            merged.rename(columns={"metaphlan_abundance_m": "metaphlan_abundance"}, inplace=True)
        if "bracken_abundance_b" in merged.columns:
            merged.rename(columns={"bracken_abundance_b": "bracken_abundance"}, inplace=True)
        # Drop helper columns
        for c in ["taxon_m", "taxon_b"]:
            if c in merged.columns:
                merged.drop(columns=[c], inplace=True)
    else:
        merged = pd.merge(meta_df, br_df, on=["sample", "taxon", "rank"], how="outer")
        # If both inputs carried a taxid column, pandas will suffix them (_x/_y).
        # Coalesce into a single 'taxid' column for downstream consumers (e.g., genome-length correction).
        if ("taxid_x" in merged.columns) or ("taxid_y" in merged.columns):
            # Prefer MetaPhlAn taxid when present, else Bracken's
            merged["taxid"] = merged.get("taxid_x")
            if "taxid_y" in merged.columns:
                merged["taxid"] = merged["taxid"].where(merged["taxid"].notna(), merged["taxid_y"]) 
            # Drop helper columns to avoid confusion
            for c in ("taxid_x", "taxid_y"):
                if c in merged.columns:
                    merged.drop(columns=[c], inplace=True)
    # Ensure numeric dtype and fill NaNs to avoid FutureWarning
    merged["metaphlan_abundance"] = pd.to_numeric(merged["metaphlan_abundance"], errors="coerce").fillna(0.0).astype(float)
    merged["bracken_abundance"] = pd.to_numeric(merged["bracken_abundance"], errors="coerce").fillna(0.0).astype(float)
    # Defensive: no negative abundances
    merged.loc[merged["metaphlan_abundance"] < 0, "metaphlan_abundance"] = 0.0
    merged.loc[merged["bracken_abundance"] < 0, "bracken_abundance"] = 0.0

    # Wide view: two columns side-by-side
    wide = merged.rename(columns={"metaphlan_abundance": "metaphlan", "bracken_abundance": "bracken"})
    wide = wide[[c for c in ["sample", "taxon", "rank", "metaphlan", "bracken"] if c in wide.columns]]

    # Sort consistently using ordered categorical rank
    if not merged.empty:
        cat = pd.Categorical(merged["rank"], categories=RANK_ORDER, ordered=True)
        merged = merged.assign(_rk=cat).sort_values(["_rk", "taxon", "sample"]).drop(columns=["_rk"]).reset_index(drop=True)
    if not wide.empty:
        catw = pd.Categorical(wide["rank"], categories=RANK_ORDER, ordered=True)
        wide = wide.assign(_rk=catw).sort_values(["_rk", "taxon", "sample"]).drop(columns=["_rk"]).reset_index(drop=True)

    return merged, wide

def validate_args(args):
    if not os.path.isdir(args.metaphlan_dir):
        raise ValueError(f"MetaPhlAn directory not found: {args.metaphlan_dir}")
    if not os.path.isdir(args.bracken_dir):
        raise ValueError(f"Bracken directory not found: {args.bracken_dir}")
    if args.pseudocount <= 0:
        raise ValueError(f"Pseudocount must be positive: {args.pseudocount}")
    # Validate normalizations list
    valid_norms = {"clr","tmm","rle","css","rcss","uq","median","total","quantile"}
    for m in args.normalizations:
        if m not in valid_norms:
            raise ValueError(f"Unknown normalization '{m}'. Valid: {sorted(valid_norms)}")
    # Guard degenerate trims
    if not (0.0 <= args.tmm_lfc_trim < 1.0) or not (0.0 <= args.tmm_a_trim < 1.0):
        raise ValueError("TMM trim parameters must be within [0,1).")
    if args.css_rel_tol < 0:
        raise ValueError("--css-rel-tol must be non-negative")

def main():
    ap = argparse.ArgumentParser(description="Merge MetaPhlAn and Bracken outputs into unified TSVs")
    ap.add_argument("--metaphlan-dir", default="metaphlan_out", help="Directory with *_profile.txt files")
    ap.add_argument("--bracken-dir", default="bracken_output", help="Directory with *.{rank}.bracken files")
    ap.add_argument("--rank", default="species", choices=list(RANK_PREFIX.keys()), help="Taxonomic rank to merge (ignored if --all-ranks is set)")
    ap.add_argument("--all-ranks", action="store_true", help="Merge all standard ranks (phylum,class,order,family,genus,species) and concatenate results")
    ap.add_argument("--sample-filter", default=None, help="Optional regex to filter sample filenames (applied to basenames)")
    ap.add_argument("--out-dir", default="diversity_results", help="Directory to place outputs (created if missing)")
    ap.add_argument("--out-long", default="metaphlan_bracken_merged.tsv", help="Output file name for tidy (long) TSV (will be placed under --out-dir unless absolute path)")
    ap.add_argument("--out-wide", default="metaphlan_bracken_merged_wide.tsv", help="Output file name for wide TSV (will be placed under --out-dir unless absolute path)")
    ap.add_argument("--out-both", default="metaphlan_bracken_both.tsv", help="Output file name for taxa present in both tools (wide format; placed under --out-dir unless absolute path)")
    ap.add_argument("--limit", type=int, default=None, help="Optional limit on number of samples processed (per tool)")
    ap.add_argument("--join-on", choices=["name", "taxid"], default="name", help="Merge key for MetaPhlAn vs Bracken: by taxon name or NCBI taxonomy id")
    ap.add_argument("--compute-clr", action="store_true", default=True, help="Compute CLR (centered log-ratio) per sample+rank for both tools using pseudocount")
    ap.add_argument("--no-clr", dest="compute_clr", action="store_false", help="Disable CLR computation")
    ap.add_argument("--pseudocount", type=float, default=1e-6, help="Zero-replacement pseudocount in proportion units; applied only to zeros (default 1e-6)")
    ap.add_argument("--bacteria-only", action="store_true", help="Keep only taxa that belong to Bacteria (based on MetaPhlAn lineages)")
    ap.add_argument("--genome-lengths", default=None, help="Optional TSV/CSV with columns: taxid, genome_length[, ploidy]. Used to approximate Bracken taxonomic abundance via genome-length correction.")
    # Normalization methods based on Pereira et al. 2018 (BMC Genomics) comparison
    ap.add_argument("--normalizations", type=lambda s: [x.strip().lower() for x in s.split(",") if x.strip()], default=["tmm","rle","css","uq","total"], help="Comma-separated list of normalization methods to compute per tool: clr,tmm,rle,css,rcss,uq,median,total,quantile. Default: tmm,rle,css,uq,total")
    ap.add_argument("--css-rel-tol", dest="css_rel_tol", type=float, default=0.1, help="Relative tolerance for CSS quantile agreement (fraction of median at that quantile); default 0.1")
    ap.add_argument("--tmm-lfc-trim", dest="tmm_lfc_trim", type=float, default=0.3, help="Proportion to trim from log-fold-change tails in TMM (e.g., 0.3 => 15% each tail)")
    ap.add_argument("--tmm-a-trim", dest="tmm_a_trim", type=float, default=0.05, help="Proportion to trim from A-values tails in TMM (e.g., 0.05 => 2.5% each tail)")
    args = ap.parse_args()

    # Basic validation
    try:
        validate_args(args)
    except Exception as e:
        print(f"ERROR: {e}")
        raise

    metaphlan_dir = args.metaphlan_dir
    bracken_dir = args.bracken_dir
    rank = args.rank
    sample_filter = args.sample_filter
    out_dir = args.out_dir

    # Ensure output directory exists
    try:
        os.makedirs(out_dir, exist_ok=True)
    except Exception as e:
        print(f"ERROR: cannot create output directory {out_dir}: {e}")
        raise

    # If limiting, wrap the globbing by trimming after merge; simpler to pass regex filter
    if args.all_ranks:
        ranks = ["phylum", "class", "order", "family", "genus", "species"]
        merged_frames = []
        for rk in ranks:
            m, _w = merge_metaphlan_bracken(metaphlan_dir, bracken_dir, rk, sample_glob=sample_filter, join_on=args.join_on)
            merged_frames.append(m)
        merged = pd.concat(merged_frames, ignore_index=True) if merged_frames else pd.DataFrame(columns=["sample","taxon","rank","metaphlan_abundance","bracken_abundance"]) 
    else:
        merged, _wide_ignored = merge_metaphlan_bracken(metaphlan_dir, bracken_dir, rank, sample_glob=sample_filter, join_on=args.join_on)

    # If limit set, reduce to first N samples present to support quick dry runs
    if args.limit is not None and args.limit > 0:
        keep_samples = list(dict.fromkeys(merged["sample"].tolist()))[: args.limit]
        merged = merged[merged["sample"].isin(keep_samples)]

    # Optional filter: keep only bacterial taxa, determined from MetaPhlAn lineages
    if args.bacteria_only and not merged.empty:
        # Build list of MetaPhlAn profile files consistent with sample_filter (text and biom)
        meta_files = sorted(set(glob.glob(os.path.join(metaphlan_dir, "*_profile.txt")) + glob.glob(os.path.join(metaphlan_dir, "*.biom"))))
        if sample_filter:
            meta_files = [p for p in meta_files if re.search(sample_filter, os.path.basename(p))]
        filter_ranks = sorted(merged["rank"].unique().tolist())
        allowed = collect_bacterial_taxa(meta_files, filter_ranks)
        mask = pd.Series(False, index=merged.index)
        for rk, taxa_set in allowed.items():
            if taxa_set:
                mask = mask | ((merged["rank"] == rk) & (merged["taxon"].isin(list(taxa_set))))
        merged = merged[mask].reset_index(drop=True)

    # Compute CLR columns if enabled
    if args.compute_clr and not merged.empty:
        # Convert to proportions; replace zeros only with pseudocount
        def compute_clr(colname: str, outname: str):
            props = (merged[colname].astype(float) / 100.0)
            props = props.mask(props <= 0.0, args.pseudocount)
            logs = np.log(props)
            mean_logs = logs.groupby([merged["sample"], merged["rank"]]).transform("mean")
            merged[outname] = (logs - mean_logs).astype(float)

        compute_clr("metaphlan_abundance", "metaphlan_clr")
        compute_clr("bracken_abundance", "bracken_clr")
        merged["delta_clr"] = (merged["bracken_clr"] - merged["metaphlan_clr"]).astype(float)

    # Compute robust CLR (rCLR) on non-zero values only (DEICODE-style), zeros -> NaN
    def compute_rclr(colname: str, outname: str):
        if merged.empty:
            return
        vals = merged[colname].astype(float) / 100.0
        out = pd.Series(np.nan, index=merged.index, dtype=float)
        grp = merged.groupby([merged["sample"], merged["rank"]])
        for _, idx in grp.groups.items():
            ix = list(idx)
            v = vals.iloc[ix]
            mask = v > 0
            if mask.sum() == 0:
                continue
            logs = np.log(v[mask])
            mean_logs = logs.mean()
            out.iloc[np.array(ix)[mask.values]] = (logs - mean_logs).astype(float).values
        merged[outname] = out

    compute_rclr("metaphlan_abundance", "metaphlan_rclr")
    compute_rclr("bracken_abundance", "bracken_rclr")

    # Optional: Approximate Bracken taxonomic abundance via genome-length correction
    if args.genome_lengths and not merged.empty:
        gl_path = args.genome_lengths
        try:
            # Read TSV/CSV flexibly
            sep = "," if gl_path.endswith((".csv", ".CSV")) else "\t"
            gl = pd.read_csv(gl_path, sep=sep)
            # Normalize columns
            cols = {c.lower(): c for c in gl.columns}
            if "taxid" not in cols or "genome_length" not in cols:
                print(f"WARN: genome-lengths file {gl_path} missing required columns 'taxid' and 'genome_length'")
            else:
                gl2 = gl.rename(columns={cols["taxid"]: "taxid", cols["genome_length"]: "genome_length"})
                if "ploidy" in cols:
                    gl2 = gl2.rename(columns={cols["ploidy"]: "ploidy"})
                else:
                    gl2["ploidy"] = 1.0
                # Coerce types
                gl2["taxid"] = gl2["taxid"].astype(str)
                gl2["genome_length"] = pd.to_numeric(gl2["genome_length"], errors="coerce")
                gl2["ploidy"] = pd.to_numeric(gl2["ploidy"], errors="coerce").fillna(1.0)
                # Attach genome lengths to merged via taxid
                merged["taxid"] = merged["taxid"].astype(str)
                merged = merged.merge(gl2[["taxid","genome_length","ploidy"]], on="taxid", how="left")
                # Compute corrected abundance per sample+rank vectorized, ignoring rows without genome_length
                ba = pd.to_numeric(merged["bracken_abundance"], errors="coerce") / 100.0
                glv = pd.to_numeric(merged["genome_length"], errors="coerce")
                plv = pd.to_numeric(merged["ploidy"], errors="coerce").fillna(1.0)
                valid = (ba > 0) & glv.notna() & (glv > 0)
                adj = pd.Series(np.nan, index=merged.index, dtype=float)
                adj.loc[valid] = (ba.loc[valid] / (glv.loc[valid] * plv.loc[valid])).astype(float)
                grp_sum = adj.groupby([merged["sample"], merged["rank"]]).transform("sum")
                with np.errstate(invalid="ignore", divide="ignore"):
                    merged["bracken_taxonomic_abundance"] = (adj / grp_sum) * 100.0
                print("NOTE: Computed bracken_taxonomic_abundance using genome-length correction (approximate taxonomic abundance)")
        except Exception as e:
            print(f"WARN: failed to apply genome-length correction: {e}")

    # Additional normalization methods per tool (TMM, RLE, CSS, RCSS, UQ, Median, Total, Quantile)
    def _pivot_rank(df: pd.DataFrame, value_col: str, rank_value: str) -> pd.DataFrame:
        sub = df[df["rank"] == rank_value][["taxon","sample", value_col]].copy()
        # Aggregate duplicates by summing abundances for same sample+taxon
        agg = sub.groupby(["taxon","sample"], as_index=False)[value_col].sum()
        mat = agg.pivot(index="taxon", columns="sample", values=value_col).fillna(0.0)
        # Ensure float
        return mat.astype(float)

    def _renorm_cols_to_percent(mat: pd.DataFrame) -> pd.DataFrame:
        col_sums = mat.sum(axis=0)
        col_sums = col_sums.replace(0.0, 1.0)
        return mat.divide(col_sums, axis=1) * 100.0

    def _add_eps(mat: pd.DataFrame, eps: float) -> pd.DataFrame:
        return mat + eps

    def norm_total(mat: pd.DataFrame) -> pd.DataFrame:
        # Identity under re-normalization
        return _renorm_cols_to_percent(mat)

    def norm_uq(mat: pd.DataFrame) -> pd.DataFrame:
        uq = mat.quantile(0.75, axis=0)
        uq = uq.replace(0.0, 1.0)
        scaled = mat.divide(uq, axis=1)
        return _renorm_cols_to_percent(scaled)

    def norm_median(mat: pd.DataFrame) -> pd.DataFrame:
        med = mat.median(axis=0)
        med = med.replace(0.0, 1.0)
        scaled = mat.divide(med, axis=1)
        return _renorm_cols_to_percent(scaled)

    def norm_rle(mat: pd.DataFrame, eps: float) -> pd.DataFrame:
        m = _add_eps(mat, eps)
        # Geometric mean across samples for each taxon
        gm = np.exp(np.log(m).mean(axis=1))
        # Avoid divide by zero: keep taxa with gm>0
        gm = gm.replace(0.0, np.nan)
        ratio = m.divide(gm, axis=0)
        factors = ratio.median(axis=0, skipna=True)
        factors = factors.replace(0.0, 1.0).fillna(1.0)
        scaled = mat.divide(factors, axis=1)
        return _renorm_cols_to_percent(scaled)

    def _tmm_factor(col: pd.Series, ref: pd.Series, lfc_trim: float, a_trim: float, eps: float) -> float:
        x = col.astype(float).values
        r = ref.astype(float).values
        x = x + eps
        r = r + eps
        M = np.log2(x) - np.log2(r)
        A = 0.5 * (np.log2(x) + np.log2(r))
        # Trim extremes
        keep = np.isfinite(M) & np.isfinite(A)
        M2 = M[keep]
        A2 = A[keep]
        if M2.size == 0:
            return 1.0
        # Determine trim bounds
        lt = lfc_trim / 2.0
        at = a_trim / 2.0
        lowM, highM = np.quantile(M2, [lt, 1.0 - lt]) if M2.size > 1 else (M2[0], M2[0])
        lowA, highA = np.quantile(A2, [at, 1.0 - at]) if A2.size > 1 else (A2[0], A2[0])
        mask = (M >= lowM) & (M <= highM) & (A >= lowA) & (A <= highA) & keep
        if not np.any(mask):
            return 1.0
        meanM = float(np.mean(M[mask]))
        return float(2.0 ** meanM)

    def norm_tmm(mat: pd.DataFrame, lfc_trim: float, a_trim: float, eps: float) -> pd.DataFrame:
        # Choose reference as column with median 75th percentile
        uq = mat.quantile(0.75, axis=0)
        ref_name = uq.sort_values().index[len(uq) // 2]
        ref = mat[ref_name]
        factors = {}
        for c in mat.columns:
            factors[c] = _tmm_factor(mat[c], ref, lfc_trim, a_trim, eps)
        factors = pd.Series(factors)
        factors = factors.replace(0.0, 1.0).fillna(1.0)
        scaled = mat.divide(factors, axis=1)
        return _renorm_cols_to_percent(scaled)

    def norm_css(mat: pd.DataFrame, rel_tol: float) -> pd.DataFrame:
        # Determine a common quantile threshold q in [0.5..0.9] where sample quantiles agree
        qs = [0.5, 0.6, 0.7, 0.8, 0.9]
        chosen_q = 0.5
        for q in qs:
            qvals = mat.quantile(q, axis=0)
            med = float(np.median(qvals.values))
            if med <= 0:
                continue
            dev = np.abs(qvals - med)
            max_dev = float(dev.max())
            if max_dev <= rel_tol * med:
                chosen_q = q
        # Compute thresholds per sample at chosen_q
        q_per_sample = mat.quantile(chosen_q, axis=0)
        # Sum of low-abundant taxa up to threshold
        factors = {}
        for c in mat.columns:
            thr = q_per_sample[c]
            col = mat[c]
            factors[c] = float(col[col <= thr].sum())
        factors = pd.Series(factors)
        factors = factors.replace(0.0, 1.0)
        # Scale so that median factor maps to 100
        median_factor = float(np.median(factors.values)) or 1.0
        scaled = mat.divide(factors, axis=1) * median_factor
        return _renorm_cols_to_percent(scaled)

    def norm_rcss(mat: pd.DataFrame) -> pd.DataFrame:
        factors = {}
        for c in mat.columns:
            thr = mat[c].median()
            col = mat[c]
            factors[c] = float(col[col >= thr].sum())
        factors = pd.Series(factors)
        factors = factors.replace(0.0, 1.0)
        median_factor = float(np.median(factors.values)) or 1.0
        scaled = mat.divide(factors, axis=1) * median_factor
        return _renorm_cols_to_percent(scaled)

    def norm_quantile(mat: pd.DataFrame) -> pd.DataFrame:
        # Standard quantile normalization across columns
        sorted_mat = np.sort(mat.values, axis=0)
        mean_sorted = np.mean(sorted_mat, axis=1)
        # Map back preserving ranks per column
        result = mat.copy()
        for j, c in enumerate(mat.columns):
            order = np.argsort(mat[c].values)
            inv_order = np.empty_like(order)
            inv_order[order] = np.arange(order.size)
            result[c] = mean_sorted[inv_order]
        return _renorm_cols_to_percent(result)

    def apply_norms(df: pd.DataFrame, methods: List[str]) -> pd.DataFrame:
        if df.empty:
            return df
        methods = [m for m in methods if m in {"tmm","rle","css","rcss","uq","median","total","quantile"}]
        if not methods:
            return df
        pct_eps = args.pseudocount * 100.0
        for rk in sorted(df["rank"].unique().tolist()):
            for tool_col, out_prefix in [("metaphlan_abundance","metaphlan"),("bracken_abundance","bracken")]:
                mat = _pivot_rank(df, tool_col, rk)
                if mat.shape[1] == 0:
                    continue
                results: Dict[str, pd.DataFrame] = {}
                for m in methods:
                    try:
                        if m == "total":
                            res = norm_total(mat)
                        elif m == "uq":
                            res = norm_uq(mat)
                        elif m == "median":
                            res = norm_median(mat)
                        elif m == "rle":
                            res = norm_rle(mat, pct_eps)
                        elif m == "tmm":
                            res = norm_tmm(mat, args.tmm_lfc_trim, args.tmm_a_trim, pct_eps)
                        elif m == "css":
                            res = norm_css(mat, args.css_rel_tol)
                        elif m == "rcss":
                            res = norm_rcss(mat)
                        elif m == "quantile":
                            res = norm_quantile(mat)
                        else:
                            continue
                        r_results = res.stack().rename(f"{out_prefix}_{m}").reset_index().rename(columns={"level_0":"taxon","level_1":"sample"})
                        # Merge back into df for this rank only
                        mask = (df["rank"] == rk)
                        df.loc[mask, f"{out_prefix}_{m}"] = pd.merge(
                            df.loc[mask, ["sample","taxon"]],
                            r_results,
                            on=["sample","taxon"],
                            how="left",
                        )[f"{out_prefix}_{m}"]
                        df[f"{out_prefix}_{m}"] = pd.to_numeric(df[f"{out_prefix}_{m}"], errors="coerce")
                    except Exception as e:
                        print(f"WARN: failed to compute {m} normalization for rank {rk} on {out_prefix}: {e}")
        return df

    # Apply requested normalizations (excluding CLR which is handled above)
    non_clr_methods = [m for m in args.normalizations if m != "clr"]
    if non_clr_methods:
        merged = apply_norms(merged, non_clr_methods)

    # Build wide view including CLR/rCLR and requested normalization columns if present
    wide = merged.rename(columns={"metaphlan_abundance": "metaphlan", "bracken_abundance": "bracken"})
    wide_cols = ["sample", "taxon", "rank", "metaphlan", "bracken"]
    # Add CLR columns if present
    if args.compute_clr and "metaphlan_clr" in merged.columns:
        wide_cols.extend(["metaphlan_clr", "bracken_clr", "delta_clr"])
    # Add rCLR columns if present
    for c in ["metaphlan_rclr","bracken_rclr"]:
        if c in merged.columns:
            wide_cols.append(c)
    # Add abundance-type clarity columns
    wide["metaphlan_taxonomic"] = wide.get("metaphlan", np.nan)
    if "bracken_taxonomic_abundance" in merged.columns:
        # Direct assignment preserves row count; avoid merging that could duplicate rows
        wide["bracken_taxonomic"] = merged["bracken_taxonomic_abundance"].values
        wide_cols.append("bracken_taxonomic")
    # Add extra normalization columns if present
    for m in non_clr_methods:
        m_meta = f"metaphlan_{m}"
        m_br = f"bracken_{m}"
        if m_meta in merged.columns and m_br in merged.columns:
            wide_cols.extend([m_meta, m_br])
    wide = wide[wide_cols]

    # Resolve output paths under out_dir unless absolute paths provided
    out_long_path = args.out_long if os.path.isabs(args.out_long) else os.path.join(out_dir, args.out_long)
    out_wide_path = args.out_wide if os.path.isabs(args.out_wide) else os.path.join(out_dir, args.out_wide)
    out_both_path = args.out_both if os.path.isabs(args.out_both) else os.path.join(out_dir, args.out_both)

    # Write outputs with high precision to avoid rounding to zero on small values
    float_format = "%.15g"  # higher precision to reduce rounding to zero
    merged.to_csv(out_long_path, sep="\t", index=False, float_format=float_format)
    wide.to_csv(out_wide_path, sep="\t", index=False, float_format=float_format)

    # Also emit a table including only taxa observed by both tools
    both_mask = (wide["metaphlan"].astype(float) > 0.0) & (wide["bracken"].astype(float) > 0.0)
    both = wide.loc[both_mask].reset_index(drop=True)
    both.to_csv(out_both_path, sep="\t", index=False, float_format=float_format)

    print(f"Wrote: {out_long_path} (rows={len(merged)})")
    print(f"Wrote: {out_wide_path} (rows={len(wide)})")
    print(f"Wrote: {out_both_path} (rows={len(both)})")
    print("NOTE: MetaPhlAn reports taxonomic abundance (marker-based); Bracken reports sequence abundance (reads).")
    if args.genome_lengths:
        print("NOTE: When genome lengths were provided, an approximate Bracken taxonomic abundance was computed via genome-length correction.")
    # UX: quick counts summary
    try:
        print(merged["rank"].value_counts().to_string())
        print(f"samples: {merged['sample'].nunique()}, taxa: {merged['taxon'].nunique()}")
    except Exception:
        pass


if __name__ == "__main__":
    main()
