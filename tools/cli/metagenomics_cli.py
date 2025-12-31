#!/usr/bin/env python3
"""
Metagenomics unified CLI

Subcommands:
- merge:    Merge MetaPhlAn and Bracken outputs to long/wide TSVs
- diagnose: Per-sample distribution diagnostics on the merged wide TSV
- ecology:  Ecology analytics (AFD, occupancy, GOF, richness, forecasts)
- all:      Run merge -> diagnose -> ecology in sequence

Config:
- Accepts a unified TOML (see config.sample.toml) with sections [merge], [diagnose], [ecology]
- CLI flags override TOML
"""
from __future__ import annotations

import argparse
import os
import sys
import shlex
import subprocess
import csv
import glob
import math
import re
import logging
from typing import Dict, List, Any, Optional, Tuple

import numpy as np
import pandas as pd


def _load_toml(path: Optional[str]) -> Dict[str, Any]:
    if not path:
        return {}
    try:
        import tomllib as _toml  # py311+
    except ModuleNotFoundError:
        try:
            import tomli as _toml  # type: ignore
        except ModuleNotFoundError:
            print("ERROR: TOML support needs Python 3.11+ (tomllib) or 'pip install tomli' on older Pythons.")
            raise SystemExit(2)
    with open(path, "rb") as f:
        return _toml.load(f)


def _kv_to_cli(d: Dict[str, Any], flag_keys: Optional[set] = None, csv_keys: Optional[set] = None) -> List[str]:
    args: List[str] = []
    flag_keys = flag_keys or set()
    csv_keys = csv_keys or set()
    for k, v in d.items():
        key = f"--{str(k).replace('_','-')}"
        if k in flag_keys:
            if bool(v):
                args.append(key)
            continue
        # Skip None or empty strings
        if v is None:
            continue
        if isinstance(v, str) and v == "":
            continue
        # Skip booleans that aren't declared as flags
        if isinstance(v, bool):
            continue
        if isinstance(v, list) and (k in csv_keys):
            if not v:
                continue
            v = ",".join(str(x) for x in v)
        args += [key, str(v)]
    return args


def _run_py(module_filename: str, argv: List[str]) -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    candidate = os.path.join(here, module_filename)
    script = candidate if os.path.exists(candidate) else module_filename
    cmd = [sys.executable or "python3", script] + argv
    print("RUN:", " ".join(shlex.quote(c) for c in cmd))
    if getattr(_run_py, "_dry_run", False):
        return 0
    return subprocess.call(cmd)

def sub_merge(args: argparse.Namespace, cfg: Dict[str, Any]) -> int:
    conf = (cfg.get("merge") or {}).copy()
    send_no_clr = False

    # TOML disables CLR
    if conf.get("compute_clr") is False:
        send_no_clr = True
        conf.pop("compute_clr", None)
    # CLI disables CLR
    if args.no_clr:
        send_no_clr = True
        conf.pop("compute_clr", None)

    # CLI overrides
    if args.metaphlan_dir: conf["metaphlan_dir"] = args.metaphlan_dir
    if args.bracken_dir: conf["bracken_dir"] = args.bracken_dir
    if args.rank: conf["rank"] = args.rank
    if args.all_ranks: conf["all_ranks"] = True
    if args.sample_filter: conf["sample_filter"] = args.sample_filter
    if args.out_dir: conf["out_dir"] = args.out_dir
    if args.out_long: conf["out_long"] = args.out_long
    if args.out_wide: conf["out_wide"] = args.out_wide
    if args.out_both: conf["out_both"] = args.out_both
    if args.limit is not None: conf["limit"] = args.limit
    if args.join_on: conf["join_on"] = args.join_on
    if args.pseudocount is not None: conf["pseudocount"] = args.pseudocount
    if args.bacteria_only: conf["bacteria_only"] = True
    if args.genome_lengths: conf["genome_lengths"] = args.genome_lengths
    if args.normalizations: conf["normalizations"] = args.normalizations
    if args.css_rel_tol is not None: conf["css_rel_tol"] = args.css_rel_tol
    if args.tmm_lfc_trim is not None: conf["tmm_lfc_trim"] = args.tmm_lfc_trim
    if args.tmm_a_trim is not None: conf["tmm_a_trim"] = args.tmm_a_trim

    # Required args sanity check
    for k in ["metaphlan_dir", "bracken_dir", "rank", "out_dir", "out_wide"]:
        if not conf.get(k):
            print(f"[merge] missing required: {k}", file=sys.stderr)
            return 2

    # ensure parent dirs exist for outputs
    for k in ("out_dir", "out_long", "out_wide", "out_both"):
        p = conf.get(k)
        if not p:
            continue
        # If it's a filename only, use out_dir; otherwise use its dirname
        d = os.path.dirname(p) if os.path.dirname(p) else conf.get("out_dir")
        if d:
            os.makedirs(d, exist_ok=True)

    # Run merge in-process to avoid external script
    try:
        return _merge_execute(conf, send_no_clr=send_no_clr)
    except Exception as e:
        print(f"[merge] ERROR: {e}", file=sys.stderr)
        return 1


def sub_diagnose(args: argparse.Namespace, cfg: Dict[str, Any]) -> int:
    conf = (cfg.get("diagnose") or {}).copy()

    # CLI overrides
    if args.input: conf["input"] = args.input
    if args.out: conf["out"] = args.out
    if args.ranks: conf["ranks"] = args.ranks
    if args.transforms: conf["transforms"] = args.transforms
    if args.alpha is not None: conf["alpha"] = args.alpha
    if args.zero_thresh is not None: conf["zero_thresh"] = args.zero_thresh

    # ensure output dir exists if provided
    out = conf.get("out")
    if out:
        d = os.path.dirname(out)
        if d:
            os.makedirs(d, exist_ok=True)

    # Run diagnostics in-process
    try:
        return _diagnose_execute(conf)
    except Exception as e:
        print(f"[diagnose] ERROR: {e}", file=sys.stderr)
        return 1



def sub_ecology(args: argparse.Namespace, cfg: Dict[str, Any]) -> int:
    conf = (cfg.get("ecology") or {}).copy()
    # CLI overrides
    mapping = {
        "merged": args.merged,
        "counts": args.counts,
        "rank": args.rank,
        "tool_pair": args.tool_pair,
        "out_dir": args.out_dir,
        "log_level": args.log_level,
        "agreement_mode": args.agreement_mode,
        "subcomp_mode": args.subcomp_mode,
        "subcomp_prevalence": args.subcomp_prevalence,
        "subcomp_topk": args.subcomp_topk,
        "zero_replacement": args.zero_replacement,
        "pseudocount": args.pseudocount,
        "zero_delta": args.zero_delta,
        "dirichlet_alpha": args.dirichlet_alpha,
        "bias_correct": args.bias_correct,
        "bias_shrink": args.bias_shrink,
        "deming_delta": args.deming_delta,
        "detect_thresh": args.detect_thresh,
        "detect_alpha": args.detect_alpha,
        "detect_fixed_t": args.detect_fixed_t,
        "target_detect": args.target_detect,
        "depths": args.depths,
        "forecast_depths": args.forecast_depths,
        "beta_grid": args.beta_grid,
        "fit_per_species_beta": args.fit_per_species_beta,
        "afd": args.afd,
        "model_selection": args.model_selection,
        "ms_folds": args.ms_folds,
        "pln_quadrature": args.pln_quadrature,
        "zinb_occupancy": args.zinb_occupancy,
        "metadata": args.metadata,
        "sample_col": args.sample_col,
        "include_regex": args.include_regex,
        "cohort_col": args.cohort_col,
        "cohort_allow": args.cohort_allow,
        "date_col": args.date_col,
        "exclude_old_before": args.exclude_old_before,
        "neg_regex": args.neg_regex,
        "n_bootstrap": args.n_bootstrap,
        "bootstrap_block_col": args.bootstrap_block_col,
        "bootstrap_max": args.bootstrap_max,
        "bootstrap_eps": args.bootstrap_eps,
        "cluster_cols": args.cluster_cols,
        "time_col": args.time_col,
        "mbb_window": args.mbb_window,
        "stratify": args.stratify,
        "batch_col": args.batch_col,
        "priority_taxa": args.priority_taxa,
        "gof_alpha": args.gof_alpha,
        "slope_min_n": args.slope_min_n,
        "mad_min_n": args.mad_min_n,
        "report_holm": args.report_holm,
        "mt_hierarchical": args.mt_hierarchical,
        "null": getattr(args, "null", None),
        "rarefy_grid": args.rarefy_grid,
        "pseudocount_grid": args.pseudocount_grid,
        "q_grid": args.q_grid,
        "seed": args.seed,
    }
    for k, v in mapping.items():
        if v is not None and v != "":
            conf[k] = v
    # Required args sanity check
    for k in ["merged", "counts", "rank", "out_dir"]:
        if not conf.get(k):
            print(f"[ecology] missing required: {k}", file=sys.stderr)
            return 2
    try:
        return _ecology_execute(conf)
    except Exception as e:
        print(f"[ecology] ERROR: {e}", file=sys.stderr)
        return 1


def sub_all(args: argparse.Namespace, cfg: Dict[str, Any]) -> int:
    # Helper: ensure all expected keys are present in Namespace
    def fill_namespace(keys, args):
        ns = vars(args).copy()
        for k in keys:
            if k not in ns:
                ns[k] = None
        return argparse.Namespace(**ns)

    # 1) merge
    m_cfg = (cfg.get("merge") or {}).copy()
    for k in ["metaphlan_dir","bracken_dir","rank","all_ranks","sample_filter","out_dir","out_long","out_wide","out_both","limit","join_on","no_clr","pseudocount","bacteria_only","genome_lengths","normalizations","css_rel_tol","tmm_lfc_trim","tmm_a_trim"]:
        if hasattr(args, k):
            v = getattr(args, k)
            if v is not None and v != "":
                m_cfg[k] = v
    out_dir = m_cfg.get("out_dir", "diversity_results")
    os.makedirs(out_dir, exist_ok=True)
    out_wide = m_cfg.get("out_wide", "metaphlan_bracken_merged_wide.tsv")
    # Avoid double out_dir prefix
    if os.path.isabs(out_wide) or os.path.dirname(out_wide):
        merged_wide_path = out_wide
    else:
        merged_wide_path = os.path.join(out_dir, out_wide)
    merge_keys = ["metaphlan_dir","bracken_dir","rank","all_ranks","sample_filter","out_dir","out_long","out_wide","out_both","limit","join_on","no_clr","pseudocount","bacteria_only","genome_lengths","normalizations","css_rel_tol","tmm_lfc_trim","tmm_a_trim"]
    rc = sub_merge(fill_namespace(merge_keys, args), {"merge": m_cfg})
    if rc != 0:
        return rc

    # 2) diagnose
    d_cfg = (cfg.get("diagnose") or {}).copy()
    d_cfg.setdefault("input", merged_wide_path)
    d_out = d_cfg.setdefault("out", os.path.join(out_dir, "diagnostics", "sample_distribution.tsv"))
    os.makedirs(os.path.dirname(d_out), exist_ok=True)
    diagnose_keys = ["input","out","ranks","transforms","alpha","zero_thresh"]
    rc = sub_diagnose(fill_namespace(diagnose_keys, args), {"diagnose": d_cfg})
    if rc != 0:
        return rc

    # 3) ecology
    e_cfg = (cfg.get("ecology") or {}).copy()
    e_cfg.setdefault("merged", merged_wide_path)
    e_dir = e_cfg.setdefault("out_dir", os.path.join(out_dir, "eco_metrics"))
    os.makedirs(e_dir, exist_ok=True)
    if not e_cfg.get("counts"):
        print("[all] You must provide --counts (or set [ecology].counts in TOML).", file=sys.stderr)
        return 2
    ecology_keys = ["merged","counts","rank","tool_pair","out_dir","log_level","agreement_mode","subcomp_mode","subcomp_prevalence","subcomp_topk","zero_replacement","pseudocount","zero_delta","dirichlet_alpha","bias_correct","bias_shrink","deming_delta","detect_thresh","detect_alpha","detect_fixed_t","target_detect","depths","forecast_depths","beta_grid","fit_per_species_beta","afd","model_selection","ms_folds","pln_quadrature","zinb_occupancy","metadata","sample_col","include_regex","cohort_col","cohort_allow","date_col","exclude_old_before","neg_regex","n_bootstrap","bootstrap_block_col","bootstrap_max","bootstrap_eps","cluster_cols","time_col","mbb_window","stratify","batch_col","priority_taxa","gof_alpha","slope_min_n","mad_min_n","report_holm","mt_hierarchical","null","rarefy_grid","pseudocount_grid","q_grid","seed"]
    rc = sub_ecology(fill_namespace(ecology_keys, args), {"ecology": e_cfg})
    return rc


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Metagenomics unified CLI")
    p.add_argument("--config", help="Unified TOML with [merge], [diagnose], [ecology] sections")
    p.add_argument("--version", action="version", version="metagenomics-cli 0.1.0")
    sp = p.add_subparsers(dest="cmd", required=True)

    # merge
    pm = sp.add_parser("merge", help="Merge MetaPhlAn and Bracken outputs")
    pm.add_argument("--metaphlan-dir")
    pm.add_argument("--bracken-dir")
    pm.add_argument("--rank", choices=["species","genus","family","order","class","phylum"], default="species")
    pm.add_argument("--all-ranks", action="store_true")
    pm.add_argument("--sample-filter")
    pm.add_argument("--out-dir")
    pm.add_argument("--out-long")
    pm.add_argument("--out-wide")
    pm.add_argument("--out-both")
    pm.add_argument("--limit", type=int)
    pm.add_argument("--join-on", choices=["name","taxid"]) 
    pm.add_argument("--no-clr", action="store_true", help="Disable CLR computation")
    pm.add_argument("--pseudocount", type=float)
    pm.add_argument("--bacteria-only", action="store_true")
    pm.add_argument("--genome-lengths")
    pm.add_argument("--normalizations", help="CSV of methods e.g. tmm,rle,css,uq,total")
    pm.add_argument("--css-rel-tol", type=float)
    pm.add_argument("--tmm-lfc-trim", type=float)
    pm.add_argument("--tmm-a-trim", type=float)

    # diagnose
    pdg = sp.add_parser("diagnose", help="Per-sample distribution diagnostics")
    pdg.add_argument("--input")
    pdg.add_argument("--out")
    pdg.add_argument("--ranks")
    pdg.add_argument("--transforms")
    pdg.add_argument("--alpha", type=float)
    pdg.add_argument("--zero-thresh", dest="zero_thresh", type=float)

    # ecology (mirror primary flags; most can be provided by TOML)
    pe = sp.add_parser("ecology", help="Ecology analytics (AFD, occupancy, GOF, richness, forecasts)")
    pe.add_argument("--merged")
    pe.add_argument("--counts")
    pe.add_argument("--rank", choices=["species","genus","family","order","class","phylum"]) 
    pe.add_argument("--tool-pair")
    pe.add_argument("--out-dir")
    pe.add_argument("--log-level")
    pe.add_argument("--agreement-mode", choices=["clr","alr","plr"]) 
    pe.add_argument("--subcomp-mode", choices=["all","prevalence","topk"])
    pe.add_argument("--subcomp-prevalence", type=float)
    pe.add_argument("--subcomp-topk", type=int)
    pe.add_argument("--zero-replacement", choices=["pseudocount","mult","dirichlet"]) 
    pe.add_argument("--pseudocount", type=float)
    pe.add_argument("--zero-delta", type=float)
    pe.add_argument("--dirichlet-alpha", type=float)
    pe.add_argument("--bias-correct", choices=["none","clr-mean","alr-deming"]) 
    pe.add_argument("--bias-shrink", type=float)
    pe.add_argument("--deming-delta", type=float)
    pe.add_argument("--detect-thresh", choices=["fixed","kFDR","neg-poisson"]) 
    pe.add_argument("--detect-alpha", type=float)
    pe.add_argument("--detect-fixed-t", type=int)
    pe.add_argument("--target-detect", type=float)
    pe.add_argument("--depths")
    pe.add_argument("--forecast-depths")
    pe.add_argument("--beta-grid")
    pe.add_argument("--fit-per-species-beta", action="store_true")
    pe.add_argument("--afd", choices=["nb","pln","both"]) 
    pe.add_argument("--model-selection", choices=["aic","cvll","stack"]) 
    pe.add_argument("--ms-folds", type=int)
    pe.add_argument("--pln-quadrature", type=int)
    pe.add_argument("--zinb-occupancy", action="store_true")
    pe.add_argument("--metadata")
    pe.add_argument("--sample-col")
    pe.add_argument("--include-regex")
    pe.add_argument("--cohort-col")
    pe.add_argument("--cohort-allow")
    pe.add_argument("--date-col")
    pe.add_argument("--exclude-old-before")
    pe.add_argument("--neg-regex")
    pe.add_argument("--n-bootstrap", type=int)
    pe.add_argument("--bootstrap-block-col")
    pe.add_argument("--bootstrap-max", type=int)
    pe.add_argument("--bootstrap-eps", type=float)
    pe.add_argument("--cluster-cols")
    pe.add_argument("--time-col")
    pe.add_argument("--mbb-window", type=int)
    pe.add_argument("--stratify", action="store_true")
    pe.add_argument("--batch-col")
    pe.add_argument("--priority-taxa")
    pe.add_argument("--gof-alpha", type=float)
    pe.add_argument("--slope-min-n", type=int)
    pe.add_argument("--mad-min-n", type=int)
    pe.add_argument("--report-holm", action="store_true")
    pe.add_argument("--mt-hierarchical", action="store_true")
    pe.add_argument("--null", choices=["none","random-multinomial","tool-permute","neutral-lite"]) 
    pe.add_argument("--rarefy-grid")
    pe.add_argument("--pseudocount-grid")
    pe.add_argument("--q-grid")
    pe.add_argument("--seed", type=int)

    # all (pipeline)
    pall = sp.add_parser("all", help="Run merge -> diagnose -> ecology")
    pall.add_argument("--metaphlan-dir")
    pall.add_argument("--bracken-dir")
    pall.add_argument("--rank", choices=["species","genus","family","order","class","phylum"]) 
    pall.add_argument("--out-dir")
    pall.add_argument("--all-ranks", action="store_true")
    pall.add_argument("--no-clr", action="store_true")
    pall.add_argument("--bacteria-only", action="store_true")
    pall.add_argument("--sample-filter")
    pall.add_argument("--counts")
    pall.add_argument("--metadata")
    return p


def main():
    parser = build_parser()
    parser.add_argument("--dry-run", action="store_true", help="Print commands, don’t execute")
    args = parser.parse_args()
    _run_py._dry_run = getattr(args, "dry_run", False)
    cfg = _load_toml(getattr(args, "config", None))

    if args.cmd == "merge":
        rc = sub_merge(args, cfg)
    elif args.cmd == "diagnose":
        rc = sub_diagnose(args, cfg)
    elif args.cmd == "ecology":
        rc = sub_ecology(args, cfg)
    else:
        rc = sub_all(args, cfg)
    raise SystemExit(rc)


# ==========================
# Embedded diversity modules
# ==========================

# ---- Merge (from merge_profiles_to_tables.py) ----

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

RANK_ORDER = ["phylum", "class", "order", "family", "genus", "species"]
BRACKEN_MIN_COLUMNS = 7
METAPHLAN_ABUNDANCE_FALLBACK_INDEX = 2

def _detect_metaphlan_header(header_line: str) -> Optional[List[str]]:
    parts = header_line.lstrip("# ").split("\t")
    required_candidates = {"clade_name", "relative_abundance"}
    if any(col in parts for col in required_candidates):
        return parts
    return None

def _canonicalize_taxon_name(name: str) -> str:
    if not name:
        return name
    s = name.replace("_", " ")
    s = re.sub(r"\s*[\(\[].*?[\)\]]\s*", " ", s)
    s = re.sub(r"^(?:Candidatus\s+|Ca\.\s+)", "", s)
    s = s.replace('"', "").strip()
    s = re.sub(r"\s+", " ", s)
    return s

def _infer_sample_id(path: str) -> str:
    base = os.path.basename(path)
    base = re.sub(r"_profile\.txt$", "", base)
    base = re.sub(r"\.(species|genus|family|order|class|phylum)\.(bracken|breport)$", "", base, flags=re.IGNORECASE)
    base = re.sub(r"\.(json|biom)$", "", base, flags=re.IGNORECASE)
    return base

def _normalize_taxon_from_metaphlan(clade_name: str, rank: str) -> Optional[str]:
    prefix = RANK_PREFIX.get(rank)
    if not prefix:
        return None
    parts = clade_name.split("|")
    for part in reversed(parts):
        if part.startswith(prefix):
            name = part[len(prefix):]
            if not name or name.lower() in {"unclassified", "unknown", "nan"}:
                return None
            return name.replace("_", " ")
    return None

def _is_json_biom(path: str) -> bool:
    try:
        with open(path, "rb") as f:
            head = f.read(8192)
            c = head.lstrip()
            return c[:1] in (b"{", b"[")
    except Exception:
        return False

def _parse_metaphlan_biom_json(profile_path: str, rank: str) -> pd.DataFrame:
    import json
    with open(profile_path, "rb") as f:
        j = json.load(f)
    sample = None
    try:
        cols = j.get("columns") or []
        if cols and isinstance(cols, list) and isinstance(cols[0], dict):
            sample = cols[0].get("id") or cols[0].get("sample_id")
    except Exception:
        sample = None
    if not sample:
        sample = _infer_sample_id(profile_path)
    rows = j.get("rows") or []
    data = j.get("data") or []
    matrix_type = j.get("matrix_type", "sparse")
    records: List[Tuple[str, str, str, float]] = []
    if matrix_type == "sparse":
        for triplet in data:
            try:
                ri, _ci, val = triplet
            except Exception:
                continue
            try:
                clade = rows[ri].get("id")
            except Exception:
                continue
            taxon = _normalize_taxon_from_metaphlan(clade, rank)
            if not taxon:
                continue
            try:
                abundance = float(val)
            except Exception:
                continue
            records.append((sample, _canonicalize_taxon_name(taxon), rank, abundance, None))
    else:
        try:
            for ri, rowvals in enumerate(data):
                if not rowvals:
                    continue
                try:
                    clade = rows[ri].get("id")
                except Exception:
                    continue
                taxon = _normalize_taxon_from_metaphlan(clade, rank)
                if not taxon:
                    continue
                try:
                    abundance = float(rowvals[0])
                except Exception:
                    continue
                records.append((sample, _canonicalize_taxon_name(taxon), rank, abundance, None))
        except Exception:
            pass
    if not records:
        return pd.DataFrame(columns=["sample","taxon","rank","metaphlan_abundance","taxid"])
    return pd.DataFrame(records, columns=["sample","taxon","rank","metaphlan_abundance","taxid"])

def _parse_metaphlan_profile(profile_path: str, rank: str) -> pd.DataFrame:
    if _is_json_biom(profile_path):
        return _parse_metaphlan_biom_json(profile_path, rank)
    sample = _infer_sample_id(profile_path)
    records: List[Tuple[str, str, str, float]] = []
    with open(profile_path, "r") as fh:
        header_cols: Optional[List[str]] = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#") and header_cols is None:
                cols = _detect_metaphlan_header(line)
                if cols is not None:
                    header_cols = cols
                continue
            cols = line.split("\t")
            if header_cols is None:
                header_cols = ["clade_name", "NCBI_tax_id", "relative_abundance"]
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
            taxon = _normalize_taxon_from_metaphlan(clade, rank)
            if not taxon:
                continue
            try:
                ra = float(ra_str)
            except ValueError:
                continue
            records.append((sample, taxon, rank, ra, taxid))
    if not records:
        return pd.DataFrame(columns=["sample","taxon","rank","metaphlan_abundance","taxid"])
    df = pd.DataFrame(records, columns=["sample","taxon","rank","metaphlan_abundance","taxid"])
    df["taxon"] = df["taxon"].map(_canonicalize_taxon_name)
    return df

def _parse_bracken_file(bracken_path: str, rank: str) -> pd.DataFrame:
    sample = _infer_sample_id(bracken_path)
    m = re.search(r"\.(species|genus|family|order|class|phylum)\.(bracken|breport)$", bracken_path, re.IGNORECASE)
    file_rank = m.group(1) if m else None
    expected_letter = BRACKEN_RANK_LETTER.get(rank)
    records: List[Tuple[str, str, str, float]] = []
    with open(bracken_path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if (row[0] == "name" and len(row) >= 7) or (row[0].lower() in {"%", "percent"} and len(row) >= 6):
                continue
            name = None; taxid = None; lvl = None; abundance = None; taxon = None
            if len(row) >= 7:
                name, taxid, lvl = row[0], row[1], row[2]
                frac_str = row[-1]
                try:
                    frac = float(frac_str)
                except ValueError:
                    continue
                abundance = frac * 100.0
                taxon = _canonicalize_taxon_name((name or "").strip())
            elif len(row) == 6:
                pct_str, _clade_reads, _taxon_reads, lvl, taxid, name = row
                try:
                    pct = float(pct_str)
                except ValueError:
                    continue
                abundance = pct
                taxon = _canonicalize_taxon_name((name or "").strip())
            else:
                continue
            if expected_letter and (lvl or "").upper() != expected_letter:
                continue
            if not name or not taxon:
                continue
            if taxon.lower() in {"unclassified", "root", "no rank"}:
                continue
            records.append((sample, taxon, rank, abundance, taxid))
    if not records:
        return pd.DataFrame(columns=["sample","taxon","rank","bracken_abundance","taxid"])
    df = pd.DataFrame(records, columns=["sample","taxon","rank","bracken_abundance","taxid"])
    if file_rank and file_rank != rank:
        print(f"WARN: {os.path.basename(bracken_path)} is '{file_rank}', but --rank is '{rank}'.")
    return df

def _collect_bacterial_taxa_from_metaphlan(profile_path: str, ranks: List[str]) -> Dict[str, set]:
    allowed: Dict[str, set] = {rk: set() for rk in ranks}
    if _is_json_biom(profile_path):
        import json
        try:
            with open(profile_path, "rb") as f:
                j = json.load(f)
            rows = j.get("rows") or []
            for row in rows:
                clade = row.get("id") or ""
                if ("k__Bacteria" not in clade) and ("d__Bacteria" not in clade):
                    continue
                for rk in ranks:
                    taxon = _normalize_taxon_from_metaphlan(clade, rk)
                    if taxon:
                        allowed[rk].add(_canonicalize_taxon_name(taxon))
        except Exception as e:
            print(f"WARN: failed to scan JSON BIOM {profile_path}: {e}")
        return allowed
    with open(profile_path, "r") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if not cols:
                continue
            clade = cols[0]
            if ("k__Bacteria" not in clade) and ("d__Bacteria" not in clade):
                continue
            for rk in ranks:
                taxon = _normalize_taxon_from_metaphlan(clade, rk)
                if taxon:
                    allowed[rk].add(_canonicalize_taxon_name(taxon))
    return allowed

def _collect_bacterial_taxa(meta_files: List[str], ranks: List[str]) -> Dict[str, set]:
    agg: Dict[str, set] = {rk: set() for rk in ranks}
    for path in meta_files:
        try:
            res = _collect_bacterial_taxa_from_metaphlan(path, ranks)
            for rk in ranks:
                agg[rk].update(res.get(rk, set()))
        except Exception as e:
            print(f"WARN: failed to scan {path} for bacterial taxa: {e}")
    return agg

def _merge_metaphlan_bracken(metaphlan_dir: str, bracken_dir: str, rank: str, sample_glob: Optional[str] = None, join_on: str = "name") -> Tuple[pd.DataFrame, pd.DataFrame]:
    meta_pattern_txt = os.path.join(metaphlan_dir, "*_profile.txt")
    meta_pattern_biom = os.path.join(metaphlan_dir, "*.biom")
    meta_pattern_json = os.path.join(metaphlan_dir, "*.json")
    meta_files = sorted(set(glob.glob(meta_pattern_txt) + glob.glob(meta_pattern_biom) + glob.glob(meta_pattern_json)))
    if sample_glob:
        meta_files = [p for p in meta_files if re.search(sample_glob, os.path.basename(p))]
    br_suffix = BRACKEN_RANK_SUFFIX.get(rank, rank)
    br_pattern1 = os.path.join(bracken_dir, f"*.{br_suffix}.bracken")
    br_pattern2 = os.path.join(bracken_dir, f"*.{br_suffix}.breport")
    all_br = sorted(set(glob.glob(br_pattern1) + glob.glob(br_pattern2)))
    if sample_glob:
        all_br = [p for p in all_br if re.search(sample_glob, os.path.basename(p))]
    chosen: Dict[str, str] = {}
    seen_types: Dict[str, set] = {}
    for p in all_br:
        s = _infer_sample_id(p)
        seen_types.setdefault(s, set()).add(os.path.splitext(p)[1].lower())
        if s not in chosen:
            chosen[s] = p
        else:
            if p.lower().endswith('.bracken') and chosen[s].lower().endswith('.breport'):
                chosen[s] = p
    for s, types in seen_types.items():
        if {'.bracken', '.breport'}.issubset(types):
            print(f"WARN: Both .bracken and .breport found for sample {s} at rank {rank}; using .bracken")
    br_files = sorted(chosen.values())
    meta_dfs: List[pd.DataFrame] = []
    for path in meta_files:
        try:
            meta_dfs.append(_parse_metaphlan_profile(path, rank))
        except Exception as e:
            print(f"WARN: failed to parse MetaPhlAn file {path}: {e}")
    br_dfs: List[pd.DataFrame] = []
    for path in br_files:
        try:
            br_dfs.append(_parse_bracken_file(path, rank))
        except Exception as e:
            print(f"WARN: failed to parse Bracken file {path}: {e}")
    meta_dfs = [df for df in meta_dfs if df is not None and not df.empty]
    br_dfs = [df for df in br_dfs if df is not None and not df.empty]
    meta_df = pd.concat(meta_dfs, ignore_index=True) if meta_dfs else pd.DataFrame(columns=["sample","taxon","rank","metaphlan_abundance","taxid"]) 
    br_df = pd.concat(br_dfs, ignore_index=True) if br_dfs else pd.DataFrame(columns=["sample","taxon","rank","bracken_abundance","taxid"]) 
    if join_on == "taxid":
        try:
            if meta_df.shape[0] > 0 and meta_df["taxid"].isna().all():
                print("WARN: MetaPhlAn profiles contain no taxids; joining on taxid may yield few or no matches.")
        except Exception:
            pass
        merged = pd.merge(
            meta_df,
            br_df,
            on=["sample", "rank", "taxid"],
            how="outer",
            suffixes=("_m", "_b"),
        )
        merged["taxon"] = merged.get("taxon_m").fillna(merged.get("taxon_b"))
        if "metaphlan_abundance_m" in merged.columns:
            merged.rename(columns={"metaphlan_abundance_m": "metaphlan_abundance"}, inplace=True)
        if "bracken_abundance_b" in merged.columns:
            merged.rename(columns={"bracken_abundance_b": "bracken_abundance"}, inplace=True)
        for c in ["taxon_m", "taxon_b"]:
            if c in merged.columns:
                merged.drop(columns=[c], inplace=True)
    else:
        merged = pd.merge(meta_df, br_df, on=["sample","taxon","rank"], how="outer")
        if ("taxid_x" in merged.columns) or ("taxid_y" in merged.columns):
            merged["taxid"] = merged.get("taxid_x")
            if "taxid_y" in merged.columns:
                merged["taxid"] = merged["taxid"].where(merged["taxid"].notna(), merged["taxid_y"]) 
            for c in ("taxid_x", "taxid_y"):
                if c in merged.columns:
                    merged.drop(columns=[c], inplace=True)
    merged["metaphlan_abundance"] = pd.to_numeric(merged["metaphlan_abundance"], errors="coerce").fillna(0.0).astype(float)
    merged["bracken_abundance"] = pd.to_numeric(merged["bracken_abundance"], errors="coerce").fillna(0.0).astype(float)
    merged.loc[merged["metaphlan_abundance"] < 0, "metaphlan_abundance"] = 0.0
    merged.loc[merged["bracken_abundance"] < 0, "bracken_abundance"] = 0.0
    wide = merged.rename(columns={"metaphlan_abundance": "metaphlan", "bracken_abundance": "bracken"})
    wide = wide[[c for c in ["sample","taxon","rank","metaphlan","bracken"] if c in wide.columns]]
    if not merged.empty:
        cat = pd.Categorical(merged["rank"], categories=RANK_ORDER, ordered=True)
        merged = merged.assign(_rk=cat).sort_values(["_rk","taxon","sample"]).drop(columns=["_rk"]).reset_index(drop=True)
    if not wide.empty:
        catw = pd.Categorical(wide["rank"], categories=RANK_ORDER, ordered=True)
        wide = wide.assign(_rk=catw).sort_values(["_rk","taxon","sample"]).drop(columns=["_rk"]).reset_index(drop=True)
    return merged, wide

def _merge_validate_args(conf: Dict[str, Any]):
    metaphlan_dir = conf.get("metaphlan_dir")
    bracken_dir = conf.get("bracken_dir")
    if not os.path.isdir(metaphlan_dir):
        raise ValueError(f"MetaPhlAn directory not found: {metaphlan_dir}")
    if not os.path.isdir(bracken_dir):
        raise ValueError(f"Bracken directory not found: {bracken_dir}")
    if conf.get("pseudocount", 1e-6) <= 0:
        raise ValueError(f"Pseudocount must be positive: {conf.get('pseudocount')}")
    if conf.get("css_rel_tol", 0.1) < 0:
        raise ValueError("--css-rel-tol must be non-negative")
    valid_norms = {"clr","tmm","rle","css","rcss","uq","median","total","quantile"}
    for m in (conf.get("normalizations") or []):
        if m not in valid_norms:
            raise ValueError(f"Unknown normalization '{m}'. Valid: {sorted(valid_norms)}")

def _merge_execute(conf: Dict[str, Any], send_no_clr: bool = False) -> int:
    conf = conf.copy()
    if send_no_clr:
        conf["compute_clr"] = False
    _merge_validate_args(conf)
    metaphlan_dir = conf["metaphlan_dir"]
    bracken_dir = conf["bracken_dir"]
    rank = conf.get("rank", "species")
    sample_filter = conf.get("sample_filter")
    out_dir = conf.get("out_dir", "diversity_results")
    os.makedirs(out_dir, exist_ok=True)
    # All-ranks option
    if conf.get("all_ranks"):
        ranks = ["phylum","class","order","family","genus","species"]
        merged_frames = []
        for rk in ranks:
            m, _ = _merge_metaphlan_bracken(metaphlan_dir, bracken_dir, rk, sample_glob=sample_filter, join_on=conf.get("join_on","name"))
            merged_frames.append(m)
        merged = pd.concat(merged_frames, ignore_index=True) if merged_frames else pd.DataFrame(columns=["sample","taxon","rank","metaphlan_abundance","bracken_abundance"]) 
    else:
        merged, _ = _merge_metaphlan_bracken(metaphlan_dir, bracken_dir, rank, sample_glob=sample_filter, join_on=conf.get("join_on","name"))
    # Limit samples if requested
    limit = conf.get("limit")
    if limit is not None and int(limit) > 0 and not merged.empty:
        keep_samples = list(dict.fromkeys(merged["sample"].tolist()))[: int(limit)]
        merged = merged[merged["sample"].isin(keep_samples)]
    # Bacteria-only filter
    if conf.get("bacteria_only") and not merged.empty:
        meta_files = sorted(set(glob.glob(os.path.join(metaphlan_dir, "*_profile.txt")) + glob.glob(os.path.join(metaphlan_dir, "*.biom"))))
        if sample_filter:
            meta_files = [p for p in meta_files if re.search(sample_filter, os.path.basename(p))]
        filter_ranks = sorted(merged["rank"].unique().tolist())
        allowed = _collect_bacterial_taxa(meta_files, filter_ranks)
        mask = pd.Series(False, index=merged.index)
        for rk, taxa_set in allowed.items():
            if taxa_set:
                mask = mask | ((merged["rank"] == rk) & (merged["taxon"].isin(list(taxa_set))))
        merged = merged[mask].reset_index(drop=True)
    # CLR
    if conf.get("compute_clr", True) and not merged.empty:
        pseudocount = float(conf.get("pseudocount", 1e-6))
        def compute_clr(colname: str, outname: str):
            props = (merged[colname].astype(float) / 100.0)
            props = props.mask(props <= 0.0, pseudocount)
            logs = np.log(props)
            mean_logs = logs.groupby([merged["sample"], merged["rank"]]).transform("mean")
            merged[outname] = (logs - mean_logs).astype(float)
        compute_clr("metaphlan_abundance", "metaphlan_clr")
        compute_clr("bracken_abundance", "bracken_clr")
        merged["delta_clr"] = (merged["bracken_clr"] - merged["metaphlan_clr"]).astype(float)
    # rCLR
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
    # Extra normalizations
    def _pivot_rank(df: pd.DataFrame, value_col: str, rank_value: str) -> pd.DataFrame:
        sub = df[df["rank"] == rank_value][["taxon","sample", value_col]].copy()
        agg = sub.groupby(["taxon","sample"], as_index=False)[value_col].sum()
        mat = agg.pivot(index="taxon", columns="sample", values=value_col).fillna(0.0)
        return mat.astype(float)
    def _renorm_cols_to_percent(mat: pd.DataFrame) -> pd.DataFrame:
        col_sums = mat.sum(axis=0)
        col_sums = col_sums.replace(0.0, 1.0)
        return mat.divide(col_sums, axis=1) * 100.0
    def _add_eps(mat: pd.DataFrame, eps: float) -> pd.DataFrame:
        return mat + eps
    def norm_total(mat: pd.DataFrame) -> pd.DataFrame:
        return _renorm_cols_to_percent(mat)
    def norm_uq(mat: pd.DataFrame) -> pd.DataFrame:
        uq = mat.quantile(0.75, axis=0)
        uq = uq.replace(0.0, 1.0)
        scaled = mat.divide(uq, axis=1)
        return _renorm_cols_to_percent(scaled)
    def norm_median(mat: pd.DataFrame) -> pd.DataFrame:
        med = mat.median(axis=0).replace(0.0, 1.0)
        scaled = mat.divide(med, axis=1)
        return _renorm_cols_to_percent(scaled)
    def norm_rle(mat: pd.DataFrame, eps: float) -> pd.DataFrame:
        m = _add_eps(mat, eps)
        gm = np.exp(np.log(m).mean(axis=1))
        gm = gm.replace(0.0, np.nan)
        ratio = m.divide(gm, axis=0)
        factors = ratio.median(axis=0, skipna=True).replace(0.0, 1.0).fillna(1.0)
        scaled = mat.divide(factors, axis=1)
        return _renorm_cols_to_percent(scaled)
    def _tmm_factor(col: pd.Series, ref: pd.Series, lfc_trim: float, a_trim: float, eps: float) -> float:
        x = col.astype(float).values + eps
        r = ref.astype(float).values + eps
        M = np.log2(x) - np.log2(r)
        A = 0.5 * (np.log2(x) + np.log2(r))
        keep = np.isfinite(M) & np.isfinite(A)
        M2 = M[keep]; A2 = A[keep]
        if M2.size == 0:
            return 1.0
        lt = lfc_trim / 2.0; at = a_trim / 2.0
        lowM, highM = np.quantile(M2, [lt, 1.0 - lt]) if M2.size > 1 else (M2[0], M2[0])
        lowA, highA = np.quantile(A2, [at, 1.0 - at]) if A2.size > 1 else (A2[0], A2[0])
        mask = (M >= lowM) & (M <= highM) & (A >= lowA) & (A <= highA) & keep
        if not np.any(mask):
            return 1.0
        meanM = float(np.mean(M[mask]))
        return float(2.0 ** meanM)
    def norm_tmm(mat: pd.DataFrame, lfc_trim: float, a_trim: float, eps: float) -> pd.DataFrame:
        uq = mat.quantile(0.75, axis=0)
        ref_name = uq.sort_values().index[len(uq) // 2]
        ref = mat[ref_name]
        factors = {c: _tmm_factor(mat[c], ref, lfc_trim, a_trim, eps) for c in mat.columns}
        factors = pd.Series(factors).replace(0.0, 1.0).fillna(1.0)
        scaled = mat.divide(factors, axis=1)
        return _renorm_cols_to_percent(scaled)
    def norm_css(mat: pd.DataFrame, rel_tol: float) -> pd.DataFrame:
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
        q_per_sample = mat.quantile(chosen_q, axis=0)
        factors = {}
        for c in mat.columns:
            thr = q_per_sample[c]
            col = mat[c]
            factors[c] = float(col[col <= thr].sum())
        factors = pd.Series(factors).replace(0.0, 1.0)
        median_factor = float(np.median(factors.values)) or 1.0
        scaled = mat.divide(factors, axis=1) * median_factor
        return _renorm_cols_to_percent(scaled)
    def norm_rcss(mat: pd.DataFrame) -> pd.DataFrame:
        factors = {}
        for c in mat.columns:
            thr = mat[c].median()
            col = mat[c]
            factors[c] = float(col[col >= thr].sum())
        factors = pd.Series(factors).replace(0.0, 1.0)
        median_factor = float(np.median(factors.values)) or 1.0
        scaled = mat.divide(factors, axis=1) * median_factor
        return _renorm_cols_to_percent(scaled)
    def norm_quantile(mat: pd.DataFrame) -> pd.DataFrame:
        sorted_mat = np.sort(mat.values, axis=0)
        mean_sorted = np.mean(sorted_mat, axis=1)
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
        pct_eps = float(conf.get("pseudocount", 1e-6)) * 100.0
        for rk in sorted(df["rank"].unique().tolist()):
            for tool_col, out_prefix in [("metaphlan_abundance","metaphlan"),("bracken_abundance","bracken")]:
                mat = _pivot_rank(df, tool_col, rk)
                if mat.shape[1] == 0:
                    continue
                for m in methods:
                    try:
                        if m == "total": res = norm_total(mat)
                        elif m == "uq": res = norm_uq(mat)
                        elif m == "median": res = norm_median(mat)
                        elif m == "rle": res = norm_rle(mat, pct_eps)
                        elif m == "tmm": res = norm_tmm(mat, float(conf.get("tmm_lfc_trim",0.3)), float(conf.get("tmm_a_trim",0.05)), pct_eps)
                        elif m == "css": res = norm_css(mat, float(conf.get("css_rel_tol", 0.1)))
                        elif m == "rcss": res = norm_rcss(mat)
                        elif m == "quantile": res = norm_quantile(mat)
                        else: continue
                        r_results = res.stack().rename(f"{out_prefix}_{m}").reset_index().rename(columns={"level_0":"taxon","level_1":"sample"})
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
    non_clr_methods = [m for m in (conf.get("normalizations") or []) if m != "clr"]
    if non_clr_methods:
        merged = apply_norms(merged, non_clr_methods)
    wide = merged.rename(columns={"metaphlan_abundance": "metaphlan", "bracken_abundance": "bracken"})
    wide_cols = ["sample","taxon","rank","metaphlan","bracken"]
    if conf.get("compute_clr", True) and "metaphlan_clr" in merged.columns:
        wide_cols.extend(["metaphlan_clr","bracken_clr","delta_clr"])
    for c in ["metaphlan_rclr","bracken_rclr"]:
        if c in merged.columns:
            wide_cols.append(c)
    wide["metaphlan_taxonomic"] = wide.get("metaphlan", np.nan)
    if "bracken_taxonomic_abundance" in merged.columns:
        wide["bracken_taxonomic"] = merged["bracken_taxonomic_abundance"].values
        wide_cols.append("bracken_taxonomic")
    for m in non_clr_methods:
        m_meta = f"metaphlan_{m}"; m_br = f"bracken_{m}"
        if m_meta in merged.columns and m_br in merged.columns:
            wide_cols.extend([m_meta, m_br])
    wide = wide[wide_cols]
    out_long = conf.get("out_long", "metaphlan_bracken_merged_long.tsv")
    out_wide = conf.get("out_wide", "metaphlan_bracken_merged_wide.tsv")
    out_both = conf.get("out_both", "metaphlan_bracken_merged_both.tsv")
    out_dir = conf.get("out_dir", "diversity_results")
    out_long_path = out_long if os.path.isabs(out_long) else os.path.join(out_dir, out_long)
    out_wide_path = out_wide if os.path.isabs(out_wide) else os.path.join(out_dir, out_wide)
    out_both_path = out_both if os.path.isabs(out_both) else os.path.join(out_dir, out_both)
    float_format = "%.15g"
    merged.to_csv(out_long_path, sep="\t", index=False, float_format=float_format)
    wide.to_csv(out_wide_path, sep="\t", index=False, float_format=float_format)
    both_mask = (wide["metaphlan"].astype(float) > 0.0) & (wide["bracken"].astype(float) > 0.0)
    both = wide.loc[both_mask].reset_index(drop=True)
    both.to_csv(out_both_path, sep="\t", index=False, float_format=float_format)
    print(f"Wrote: {out_long_path} (rows={len(merged)})")
    print(f"Wrote: {out_wide_path} (rows={len(wide)})")
    print(f"Wrote: {out_both_path} (rows={len(both)})")
    print("NOTE: MetaPhlAn reports taxonomic abundance (marker-based); Bracken reports sequence abundance (reads).")
    try:
        print(merged["rank"].value_counts().to_string())
        print(f"samples: {merged['sample'].nunique()}, taxa: {merged['taxon'].nunique()}")
    except Exception:
        pass
    return 0

# ---- Diagnostics (from per_sample_distribution.py) ----

def _skew_kurtosis(x: np.ndarray) -> Tuple[float, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 3:
        return float("nan"), float("nan")
    mu = float(np.mean(x))
    v = x - mu
    m2 = float(np.mean(v ** 2))
    if m2 <= 0:
        return float("nan"), float("nan")
    m3 = float(np.mean(v ** 3))
    m4 = float(np.mean(v ** 4))
    g1 = m3 / (m2 ** 1.5)
    g2 = m4 / (m2 ** 2) - 3.0
    return g1, g2

def _jarque_bera(x: np.ndarray) -> Tuple[float, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 8:
        return float("nan"), float("nan")
    g1, g2 = _skew_kurtosis(x)
    if not np.isfinite(g1) or not np.isfinite(g2):
        return float("nan"), float("nan")
    jb = n / 6.0 * (g1 ** 2 + (g2 ** 2) / 4.0)
    p = math.exp(-jb / 2.0)
    return float(jb), float(p)

def _detect_tool_columns(df: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    mapping: Dict[str, Dict[str, str]] = {}
    tools = []
    if "metaphlan" in df.columns:
        tools.append("metaphlan")
    if "bracken" in df.columns:
        tools.append("bracken")
    if "bracken_taxonomic" in df.columns:
        tools.append("bracken_taxonomic")
    for t in tools:
        mapping[t] = {}
        mapping[t]["raw"] = t
        clr_col = f"{t}_clr" if t != "bracken_taxonomic" else None
        rclr_col = f"{t}_rclr" if t != "bracken_taxonomic" else None
        if clr_col and clr_col in df.columns:
            mapping[t]["clr"] = clr_col
        if rclr_col and rclr_col in df.columns:
            mapping[t]["rclr"] = rclr_col
        mapping[t]["logp"] = t
    return mapping

def _build_summary(df: pd.DataFrame, rank: str, tool: str, col: str, transform: str, alpha: float, zero_thresh: float) -> pd.DataFrame:
    sub = df[df["rank"] == rank]
    if sub.empty or col not in sub.columns:
        return pd.DataFrame(columns=["sample","rank","tool","transform","n_taxa","n_nonzero","zero_fraction","skewness","excess_kurtosis","jb","pvalue","classification","suggested_correlation"])
    rows = []
    for sample, g in sub.groupby("sample"):
        vals = pd.to_numeric(g[col], errors="coerce")
        vals = vals[np.isfinite(vals)]
        n_taxa = int(vals.size)
        if n_taxa == 0:
            continue
        if transform == "raw":
            x = vals.values
            zeros = np.sum(x <= 0)
        elif transform == "logp":
            p = (vals / 100.0).values
            mask = p > 0
            x = np.log10(p[mask])
            zeros = np.sum(~mask)
        elif transform in ("clr","rclr"):
            x = vals.values
            zeros = int(np.sum(~np.isfinite(x)))
            x = x[np.isfinite(x)]
        else:
            continue
        n_nonzero = int(n_taxa - zeros)
        zero_fraction = float(zeros / n_taxa)
        g1, g2 = _skew_kurtosis(x)
        jb, p = _jarque_bera(x)
        is_param = (p >= alpha) and (zero_fraction < zero_thresh)
        classification = "parametric" if is_param else "nonparametric"
        suggested = "pearson" if is_param else "spearman"
        rows.append({
            "sample": sample,
            "rank": rank,
            "tool": tool,
            "transform": transform,
            "n_taxa": n_taxa,
            "n_nonzero": n_nonzero,
            "zero_fraction": zero_fraction,
            "skewness": g1,
            "excess_kurtosis": g2,
            "jb": jb,
            "pvalue": p,
            "classification": classification,
            "suggested_correlation": suggested,
        })
    return pd.DataFrame(rows)

def _diagnose_execute(conf: Dict[str, Any]) -> int:
    inp = conf.get("input", "diversity_results/metaphlan_bracken_merged_wide.tsv")
    out = conf.get("out", "diversity_results/diagnostics/sample_distribution.tsv")
    ranks = conf.get("ranks", "all")
    transforms = conf.get("transforms", "raw,clr,rclr,logp")
    alpha = float(conf.get("alpha", 0.05))
    zero_thresh = float(conf.get("zero_thresh", 0.2))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df = pd.read_csv(inp, sep="\t")
    if df.empty:
        print(f"ERROR: input is empty: {inp}")
        return 2
    tools = _detect_tool_columns(df)
    if not tools:
        print("ERROR: no tool columns detected in input.")
        return 2
    rank_list = [r.strip() for r in (df["rank"].unique().tolist() if ranks == "all" else str(ranks).split(","))]
    transform_list = [t.strip() for t in str(transforms).split(",") if t.strip()]
    out_frames: List[pd.DataFrame] = []
    for rk in rank_list:
        for tool, txmap in tools.items():
            for tr in transform_list:
                col = txmap.get(tr)
                if not col:
                    continue
                try:
                    out_frames.append(_build_summary(df, rk, tool, col, tr, alpha, zero_thresh))
                except Exception as e:
                    print(f"WARN: failed summary for rank={rk}, tool={tool}, transform={tr}: {e}")
    out_df = pd.concat(out_frames, ignore_index=True) if out_frames else pd.DataFrame()
    if not out_df.empty:
        cols = [
            "sample","rank","tool","transform","n_taxa","n_nonzero","zero_fraction",
            "skewness","excess_kurtosis","jb","pvalue","classification","suggested_correlation"
        ]
        out_df = out_df[cols]
        out_df.to_csv(out, sep="\t", index=False)
        print(f"Wrote: {out} (rows={len(out_df)})")
    else:
        print("WARN: nothing to write.")
    return 0


# ---- Ecology (embedded from p4_meta_postmerge.py; trimmed but feature-parity for NB path) ----

# Numeric guards
_LOG_MAX = 745.0
_EPS_DIV = 1e-12
_EPS_LOG = 1e-300

def _eco_setup_logging(level: str = "INFO") -> logging.Logger:
    lvl = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(level=lvl, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    return logging.getLogger("p4_eco")

def _eco_ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def _eco_parse_csv_list(s: Optional[str]) -> Optional[List[float]]:
    if not s:
        return None
    out: List[float] = []
    for part in str(s).split(","):
        part = part.strip()
        if part:
            try:
                out.append(float(part))
            except Exception:
                pass
    return out or None

def _eco_read_table_any(path: str) -> pd.DataFrame:
    pl = path.lower()
    if pl.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t")
    if pl.endswith(".csv"):
        return pd.read_csv(path)
    if pl.endswith((".xlsx", ".xls")):
        try:
            return pd.read_excel(path)
        except Exception as e:
            raise RuntimeError(f"Reading Excel requires openpyxl/xlrd: {e}")
    raise ValueError(f"Unsupported file type: {path}")

def _eco_safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return float('nan')
    xv = x[mask]; yv = y[mask]
    sx = float(np.std(xv)); sy = float(np.std(yv))
    if sx < _EPS_DIV or sy < _EPS_DIV:
        return float('nan')
    return float(np.corrcoef(xv, yv)[0, 1])

def _eco_tag_negatives_from_metadata(meta: pd.DataFrame) -> pd.Series:
    neg_tokens = ("neg", "blank", "control", "ntc")
    neg_cols = [c for c in meta.columns if any(tok in str(c).lower() for tok in neg_tokens)]
    is_neg = pd.Series(False, index=meta.index)
    for col in meta.columns:
        vals = meta[col].astype(str).str.lower()
        is_neg = is_neg | vals.str.contains("|".join(neg_tokens), regex=True, na=False)
    for col in neg_cols:
        if meta[col].dtype == bool:
            is_neg = is_neg | meta[col]
    return is_neg

def _eco_build_sample_lists(
    merged: pd.DataFrame,
    metadata_path: Optional[str],
    sample_col: Optional[str],
    include_regex: Optional[str],
    cohort_col: Optional[str],
    cohort_allow: Optional[str],
    date_col: Optional[str],
    exclude_old_before: Optional[str],
    neg_regex: str,
    out_dir: str,
    logger: logging.Logger
) -> Tuple[List[str], List[str], List[str]]:
    samples = sorted(merged["sample"].astype(str).unique())
    included = set(samples)
    negatives: set = set()
    excluded: set = set()

    if include_regex:
        patt = re.compile(include_regex)
        keep = {s for s in samples if patt.search(s)}
        drop = set(samples) - keep
        included &= keep
        excluded |= drop

    meta = None
    if metadata_path and os.path.exists(metadata_path):
        meta = _eco_read_table_any(metadata_path)
        if not sample_col:
            candidates = [c for c in meta.columns if str(c).strip().lower() in {"sample","sample_id","id","sampleid","sample code","sample_code"}]
            if candidates:
                sample_col = candidates[0]
        if not sample_col or sample_col not in meta.columns:
            sample_col = "SAMPLE CODE" if "SAMPLE CODE" in meta.columns else None

        if sample_col:
            key_to_samples: Dict[str, List[str]] = {}
            for s in samples:
                key = s if sample_col != "SAMPLE CODE" else s.split("_")[0]
                key_to_samples.setdefault(key, []).append(s)

            if cohort_col and cohort_col in meta.columns and cohort_allow:
                allowed = {str(x).strip() for x in cohort_allow.split(",")}
                disallowed_keys = set()
                for _, r in meta.iterrows():
                    key = str(r[sample_col]).strip() if sample_col in r else None
                    if key is None:
                        continue
                    if str(r[cohort_col]).strip() not in allowed:
                        disallowed_keys.add(key)
                for k in disallowed_keys:
                    for s in key_to_samples.get(k, []):
                        if s in included:
                            included.remove(s); excluded.add(s)

            if date_col and date_col in meta.columns and exclude_old_before:
                try:
                    cutoff = pd.to_datetime(exclude_old_before)
                    dates = pd.to_datetime(meta[date_col], errors="coerce")
                    old_keys = set(meta.loc[(dates.notna()) & (dates < cutoff), sample_col].astype(str))
                    for k in old_keys:
                        for s in key_to_samples.get(k, []):
                            if s in included:
                                included.remove(s); excluded.add(s)
                except Exception as e:
                    logger.warning("Date filter ignored: %s", e)

            neg_mask = _eco_tag_negatives_from_metadata(meta)
            neg_keys = set(meta.loc[neg_mask, sample_col].astype(str)) if sample_col in meta.columns else set()
            for k in neg_keys:
                for s in key_to_samples.get(k, []):
                    negatives.add(s)
        else:
            logger.warning("Could not detect sample column in metadata; skipping cohort/date filters.")
    else:
        if metadata_path:
            logger.warning("Metadata file not found: %s", metadata_path)

    if neg_regex:
        patt = re.compile(neg_regex)
        for s in samples:
            if patt.search(s):
                negatives.add(s)

    included -= negatives

    def write_list(lst: List[str], fname: str):
        pd.DataFrame({"sample": sorted(lst)}).to_csv(os.path.join(out_dir, fname), sep="\t", index=False)

    write_list(list(included), "samples_included.tsv")
    write_list(sorted(list(negatives)), "samples_negatives.tsv")
    write_list(sorted(list(excluded)), "samples_excluded.tsv")

    return (sorted(included), sorted(negatives), sorted(excluded))

def _eco_zero_replace_pseudocount(p: np.ndarray, pseudocount: float) -> np.ndarray:
    p = np.asarray(p, float)
    p = np.where(p > 0, p, float(pseudocount))
    s = p.sum()
    return p / (s if s > 0 else 1.0)

def _eco_zero_replace_mult(p: np.ndarray, delta: float = 1e-6) -> np.ndarray:
    p = np.asarray(p, float)
    D = p.size
    k = int((p <= 0).sum())
    if k == 0:
        s = p.sum(); return p / (s if s>0 else 1.0)
    d = min(float(delta), 1.0/(10.0*max(D,1)))
    z = (p <= 0)
    nz = ~z
    mass = k * d
    out = p.copy()
    out[z] = d
    if nz.any():
        scale = (1.0 - mass) / max(p[nz].sum(), 1e-15)
        out[nz] = p[nz] * scale
    s = out.sum()
    return out / (s if s > 0 else 1.0)

def _eco_zero_replace_dirichlet(n: Optional[np.ndarray], N: Optional[float], D: int, alpha: float = 0.5) -> np.ndarray:
    if n is None or N is None or N <= 0:
        return np.full(D, 1.0 / max(D,1))
    n = np.asarray(n, float)
    denom = float(N + D * alpha)
    return (n + alpha) / (denom if denom > 0 else 1.0)

def _eco_clr_from_p(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    with np.errstate(divide="ignore"):
        logp = np.log(np.clip(p, _EPS_LOG, None))
    return logp - np.mean(logp)

def _eco_alr_from_p(p: np.ndarray, ref_idx: int) -> np.ndarray:
    p = np.asarray(p, float)
    with np.errstate(divide="ignore"):
        lp = np.log(np.clip(p, _EPS_LOG, None))
    return lp - lp[int(ref_idx)]

def _eco_plr_vector_pair_from_p(pa: np.ndarray, pb: np.ndarray, rng=None, maxpairs=0):
    a = np.asarray(pa, float)
    b = np.asarray(pb, float)
    idx = np.where((a > 0) & (b > 0))[0]
    if idx.size < 2:
        return np.array([]), np.array([])
    pairs = np.array(np.triu_indices(idx.size, k=1)).T
    if maxpairs and pairs.shape[0] > maxpairs:
        rng = np.random.default_rng(17) if rng is None else rng
        sel = rng.choice(pairs.shape[0], size=maxpairs, replace=False)
        pairs = pairs[sel]
    v1 = np.log(a[idx[pairs[:,0]]]) - np.log(a[idx[pairs[:,1]]])
    v2 = np.log(b[idx[pairs[:,0]]]) - np.log(b[idx[pairs[:,1]]])
    return v1, v2

def _eco_transform_pair(a_pct: np.ndarray, b_pct: np.ndarray, mode: str,
                        zero_replacement: str, pseudocount: float,
                        alr_ref: Optional[int], alr_topk: int, plr_maxpairs: int, rng=None,
                        a_counts: Optional[np.ndarray]=None, b_counts: Optional[np.ndarray]=None, N_a: Optional[float]=None, N_b: Optional[float]=None,
                        zero_delta: float=1e-6, dirichlet_alpha: float=0.5):
    a = np.asarray(a_pct, float) * 0.01
    b = np.asarray(b_pct, float) * 0.01
    D = a.size
    if zero_replacement == "pseudocount":
        pa = _eco_zero_replace_pseudocount(a, pseudocount)
        pb = _eco_zero_replace_pseudocount(b, pseudocount)
    elif zero_replacement == "mult":
        pa = _eco_zero_replace_mult(a, zero_delta)
        pb = _eco_zero_replace_mult(b, zero_delta)
    else:
        pa = _eco_zero_replace_dirichlet(a_counts, N_a, D, dirichlet_alpha)
        pb = _eco_zero_replace_dirichlet(b_counts, N_b, D, dirichlet_alpha)
        pa = pa / max(pa.sum(), 1e-15); pb = pb / max(pb.sum(), 1e-15)
    if mode == "clr":
        v1 = _eco_clr_from_p(pa); v2 = _eco_clr_from_p(pb)
    elif mode == "alr":
        if alr_ref is not None:
            ref = int(alr_ref)
        elif alr_topk and alr_topk > 0:
            ref = int(np.argsort(pb)[-int(alr_topk):][0])
        else:
            ref = int(np.nanargmax(pb))
        v1 = _eco_alr_from_p(pa, ref); v2 = _eco_alr_from_p(pb, ref)
    else:
        v1, v2 = _eco_plr_vector_pair_from_p(pa, pb, rng=rng, maxpairs=plr_maxpairs)
    return v1, v2

def _eco_build_fixed_subcomposition(df, samples, t1, t2, mode="prevalence", prevalence=0.5, topk=0):
    pos = df[df["sample"].isin(samples)]
    sp = pos.groupby("taxon")["sample"].nunique()
    S = len(set(samples))
    prev = sp / max(S, 1)
    if mode == "all":
        keep = pos["taxon"].unique()
    elif mode == "prevalence":
        keep = prev[prev >= float(prevalence)].index.values
    else:
        if topk <= 0:
            keep = pos["taxon"].unique()
        else:
            means = pos.groupby("taxon")[[t1, t2]].mean().mean(axis=1)
            rank = (pd.DataFrame({"prev": prev}).join(means.rename("mean")).fillna(0.0)
                    .sort_values(["prev","mean"], ascending=[False,False]))
            keep = rank.index.values[:int(topk)]
    return np.asarray(sorted(set(map(str, keep))))

def _eco_build_aligned_vectors(g, taxa_keep, tool1, tool2):
    idx_map = {t: i for i, t in enumerate(taxa_keep)}
    a = np.zeros(len(taxa_keep), float)
    b = np.zeros(len(taxa_keep), float)
    tmp = g[["taxon", tool1, tool2]].copy()
    tmp[tool1] = pd.to_numeric(tmp[tool1], errors="coerce").astype(float)
    tmp[tool2] = pd.to_numeric(tmp[tool2], errors="coerce").astype(float)
    gg = tmp.groupby("taxon", as_index=False).mean(numeric_only=True).set_index("taxon")
    t_in = [t for t in taxa_keep if t in gg.index]
    if t_in:
        ixs = [idx_map[t] for t in t_in]
        sel = gg.loc[t_in]
        a[ixs] = sel[tool1].fillna(0.0).values
        b[ixs] = sel[tool2].fillna(0.0).values
    return a, b

def _eco_estimate_clr_bias(df, samples, taxa_keep, t1, t2, zero_method, pc, zero_delta, dirichlet_alpha, lam=0.1):
    diffs = []
    for s, g in df[df["sample"].isin(samples)].groupby("sample", sort=False):
        a_pct, b_pct = _eco_build_aligned_vectors(g, taxa_keep, t1, t2)
        v1, v2 = _eco_transform_pair(a_pct, b_pct, mode="clr",
                                     zero_replacement=zero_method, pseudocount=pc,
                                     alr_ref=None, alr_topk=0, plr_maxpairs=0,
                                     zero_delta=zero_delta, dirichlet_alpha=dirichlet_alpha)
        if v1.size:
            diffs.append(v2 - v1)
    if not diffs:
        return pd.DataFrame({"taxon": taxa_keep, "bias_clr": np.zeros(len(taxa_keep)), "se": np.zeros(len(taxa_keep))})
    D = np.vstack(diffs)
    mu = np.nanmean(D, axis=0)
    se = np.nanstd(D, axis=0, ddof=1) / np.sqrt(max(D.shape[0],1))
    w = (se**2) / (se**2 + lam)
    bias = w * mu
    return pd.DataFrame({"taxon": taxa_keep, "bias_clr": bias, "se": se})

def _eco_deming_regression(x, y, delta=1.0):
    x = np.asarray(x, float); y = np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]; y = y[mask]
    if x.size < 3: return 0.0, 1.0
    xbar, ybar = x.mean(), y.mean()
    Sxx = np.mean((x - xbar)**2); Syy = np.mean((y - ybar)**2); Sxy = np.mean((x - xbar)*(y - ybar))
    if abs(Sxy) < 1e-15:
        return ybar - xbar, 1.0
    b = (Syy - delta*Sxx + np.sqrt((Syy - delta*Sxx)**2 + 4*delta*Sxy**2)) / (2*Sxy)
    a = ybar - b * xbar
    return a, b

def _eco_agreement_per_sample_fixed(df, tool1, tool2, taxa_keep, mode="clr",
                                    zero_replacement="pseudocount", pseudocount=1e-6,
                                    alr_ref=None, alr_topk=0, plr_maxpairs=0, rng=None,
                                    zero_delta=1e-6, dirichlet_alpha=0.5,
                                    bias_df: Optional[pd.DataFrame]=None, bias_mode="none", deming_delta=1.0):
    rows = []
    taxa_keep = np.asarray(taxa_keep, dtype=str)
    bias_map = {}
    if bias_df is not None and "bias_clr" in bias_df.columns:
        bias_map = dict(zip(bias_df["taxon"].astype(str), bias_df["bias_clr"].astype(float)))
    for s, g in df.groupby("sample", sort=False):
        a_pct, b_pct = _eco_build_aligned_vectors(g, taxa_keep, tool1, tool2)
        if bias_mode == "clr-mean" and mode == "clr" and bias_map:
            v1c, v2c = _eco_transform_pair(a_pct, b_pct, "clr", zero_replacement, pseudocount,
                                           alr_ref, alr_topk, plr_maxpairs, rng,
                                           zero_delta=zero_delta, dirichlet_alpha=dirichlet_alpha)
            adj = np.array([bias_map.get(t, 0.0) for t in taxa_keep], float)
            v2c = v2c - adj
            v1, v2 = v1c, v2c
        elif bias_mode == "alr-deming" and mode == "alr":
            ref = int(alr_ref) if alr_ref is not None else int(np.nanargmax(b_pct))
            v1 = _eco_alr_from_p(_eco_zero_replace_pseudocount(a_pct*0.01, 1e-9), ref)
            v2 = _eco_alr_from_p(_eco_zero_replace_pseudocount(b_pct*0.01, 1e-9), ref)
            a_int, b_slope = _eco_deming_regression(v1, v2, delta=deming_delta)
            v2 = (v2 - a_int) / max(b_slope, 1e-12)
        else:
            v1, v2 = _eco_transform_pair(a_pct, b_pct, mode, zero_replacement, pseudocount,
                                         alr_ref, alr_topk, plr_maxpairs, rng,
                                         zero_delta=zero_delta, dirichlet_alpha=dirichlet_alpha)
        if v1.size < 2 or v2.size < 2 or v1.size != v2.size:
            pear = np.nan; spear = np.nan; aitch = np.nan
        else:
            pear = _eco_safe_corr(v1, v2)
            from scipy import stats as _st
            spear = float(_st.spearmanr(v1, v2, nan_policy="omit").correlation)
            aitch = float(np.linalg.norm(v1 - v2)) if mode == "clr" else np.nan
        rows.append({"sample": s, "pearson": pear, "spearman": spear, "aitchison_clr_norm": aitch})
    return pd.DataFrame(rows).sort_values("sample", kind="mergesort").reset_index(drop=True)

def _eco_build_counts_from_pct(df: pd.DataFrame, counts_df: pd.DataFrame, tool_col: str) -> pd.DataFrame:
    m = df[["sample", "taxon", tool_col]].merge(counts_df[["sample", "N"]], on="sample", how="inner")
    rows = []
    for s, g in m.groupby("sample", sort=False):
        N = int(round(float(g["N"].iloc[0])))
        p = np.clip(pd.to_numeric(g[tool_col], errors="coerce").astype(float).values, 0.0, 100.0) * 0.01
        raw = p * N
        flo = np.floor(raw)
        give = int(N - flo.sum())
        if give > 0:
            idx = np.argsort(raw - flo)[-give:]
            flo[idx] += 1
        n = flo.astype(int)
        for (i, (_, r)) in enumerate(g.iterrows()):
            rows.append({"sample": s, "taxon": r["taxon"], "N": N, "n": int(n[i])})
    return pd.DataFrame(rows)

def _eco_unbiased_moments(counts: pd.DataFrame) -> pd.DataFrame:
    out = []
    for tax, g in counts.groupby("taxon", sort=False):
        Ns = g["N"].astype(float).values
        ns = g["n"].astype(float).values
        D1 = float(np.sum(Ns))
        D2 = float(np.sum(Ns * (Ns - 1)))
        D3 = float(np.sum(Ns * (Ns - 1) * (Ns - 2)))
        if D1 <= 0:
            continue
        mu1 = float(np.sum(ns) / D1)
        mu2 = float(np.sum(ns * (ns - 1)) / D2) if D2 > 0 else np.nan
        out.append({"taxon": tax, "mu1": mu1, "mu2": mu2})
    return pd.DataFrame(out).sort_values("taxon", kind="mergesort").reset_index(drop=True)

def _eco_nb_logpmf_vec(k, mu, beta):
    from scipy import special as _sp
    k = np.asarray(k, dtype=np.int64)
    mu_safe = np.clip(mu, _EPS_LOG, None)
    b = float(beta)
    with np.errstate(divide='ignore', invalid='ignore'):
        return (_sp.gammaln(k + b) - _sp.gammaln(b) - _sp.gammaln(k + 1) +
                b * (np.log(b) - np.log(b + mu_safe)) +
                k * (np.log(mu_safe) - np.log(b + mu_safe)))

def _eco_safe_nb_zero_prob(mu, beta):
    mu = np.asarray(mu, float)
    if np.isinf(beta):
        return np.exp(-np.clip(mu, 0.0, _LOG_MAX))
    z = -float(beta) * np.log1p(np.clip(mu / max(float(beta), 1e-12), 0.0, 1e12))
    z = np.clip(z, -_LOG_MAX, _LOG_MAX)
    p0 = np.exp(z)
    return float(p0) if p0.ndim == 0 else p0

def _eco_occupancy_from_beta(bar_x: np.ndarray, Ns: np.ndarray, beta: float) -> np.ndarray:
    mu = np.outer(bar_x, Ns)
    p0 = _eco_safe_nb_zero_prob(mu, beta)
    return np.mean(1.0 - p0, axis=1)

def _eco_fit_shared_beta(bar_x: np.ndarray, Ns: np.ndarray, occ_emp: np.ndarray, grid: Optional[List[float]] = None) -> float:
    def mse(b: float) -> float:
        b = max(float(b), 1e-6)
        pred = _eco_occupancy_from_beta(bar_x, Ns, b)
        d = pred - occ_emp
        return float(np.nanmean(d*d))
    if grid:
        vals = np.array([mse(b) for b in grid], float)
        return float(grid[int(np.argmin(vals))])
    from scipy import optimize as _opt
    res = _opt.minimize_scalar(mse, bounds=(1e-6, 1e6), method="bounded")
    return float(res.x if res.success else 1.0)

def _eco_chisq_gof_nb(counts: pd.DataFrame, bar_x: pd.Series, beta_shared: float) -> pd.DataFrame:
    from scipy import stats as _st
    rows = []
    for tax, g in counts.groupby("taxon", sort=False):
        Ns = g.set_index("sample")["N"].astype(float)
        ks = g.set_index("sample")["n"].astype(int)
        mu = bar_x.loc[tax] * Ns
        obs_k, obs_cnt = np.unique(ks.values, return_counts=True)
        exp_cnt = []
        for k in obs_k:
            pk = np.exp(np.clip(_eco_nb_logpmf_vec(np.array([k]), mu.values, beta_shared)[0], -_LOG_MAX, 0))
            exp_cnt.append(float(np.sum(pk)))
        obs = obs_cnt.astype(float).tolist()
        exp = exp_cnt
        i = 0
        while i < len(exp):
            if exp[i] >= 5 or len(exp) == 1:
                i += 1; continue
            if i + 1 < len(exp):
                exp[i + 1] += exp[i]; obs[i + 1] += obs[i]; del exp[i]; del obs[i]
            else:
                exp[i - 1] += exp[i]; obs[i - 1] += obs[i]; del exp[i]; del obs[i]; i -= 1
        df = len(exp) - 1 - 1
        if df <= 0:
            rows.append({"taxon": tax, "gof_chisq": np.nan, "gof_df": 0, "gof_p": np.nan})
            continue
        chisq = float(np.sum((np.array(obs) - np.array(exp)) ** 2 / np.maximum(np.array(exp), 1e-12)))
        p = float(_st.chi2.sf(chisq, df))
        rows.append({"taxon": tax, "gof_chisq": chisq, "gof_df": df, "gof_p": p})
    return pd.DataFrame(rows).sort_values("taxon", kind="mergesort").reset_index(drop=True)

def _eco_bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, float)
    m = np.isfinite(p).sum()
    q = np.full_like(p, np.nan, dtype=float)
    if m == 0: return q
    vals = p[np.isfinite(p)]
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    res = np.empty_like(vals); res[order] = adj
    q[np.isfinite(p)] = res
    return q

def _eco_holm_stepdown(pvals):
    p = np.asarray(pvals, float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = np.isfinite(p)
    m = int(mask.sum())
    if m == 0:
        return out
    idx = np.argsort(p[mask])
    sorted_p = p[mask][idx]
    adj = np.maximum.accumulate((m - np.arange(m)) * sorted_p)
    adj = np.minimum(adj, 1.0)
    tmp = np.empty_like(sorted_p)
    tmp[idx] = adj
    out[mask] = tmp
    return out

def _eco_n_star_q(bar_x: float, beta: float, q: float) -> float:
    q = float(np.clip(q, 1e-12, 1-1e-12))
    if math.isinf(beta):
        return -math.log(1.0 - q) / max(bar_x, 1e-15)
    return float(beta * ((1 - q) ** (-1.0 / beta) - 1.0) / max(bar_x, 1e-15))

def _eco_richness_and_slope(bar_x: np.ndarray, beta: float, N: float) -> Tuple[float,float]:
    mu = bar_x * N
    p0 = _eco_safe_nb_zero_prob(mu, beta)
    R = float(np.sum(1.0 - p0))
    if np.isinf(beta):
        Rp = float(np.sum(bar_x * np.exp(-mu)))
    else:
        Rp = float(np.sum(bar_x * (1.0 + mu / beta) ** (-(beta + 1.0))))
    return R, Rp

def _eco_rarefy_counts(counts: pd.DataFrame, depths: List[float], draws=50, rng=None) -> pd.DataFrame:
    rng = np.random.default_rng(9) if rng is None else rng
    rows = []
    for s, g in counts.groupby("sample", sort=False):
        N = int(g["N"].iloc[0]); ns = g.set_index("taxon")["n"].astype(int)
        if N <= 0: continue
        for d in depths:
            d = int(min(N, max(1, int(d))))
            for _ in range(draws):
                if d == N:
                    n2 = ns.values
                else:
                    p = pd.to_numeric(ns.values, errors="coerce").astype(float)
                    p[~np.isfinite(p)] = 0.0
                    p[p < 0] = 0.0
                    tot = float(p.sum())
                    if tot <= 0:
                        n2 = np.zeros_like(p, dtype=int)
                    else:
                        p = p / tot
                        p = np.maximum(p, 0.0)
                        p = p / max(p.sum(), 1.0)
                        n2 = rng.multinomial(d, p)
                rows.append({"sample": s, "depth": float(d), "richness": int((n2>0).sum())})
    return pd.DataFrame(rows)

def _ecology_execute(conf: Dict[str, Any]) -> int:
    merged_path = conf.get("merged")
    counts_path = conf.get("counts")
    rank = conf.get("rank", "species")
    tool_pair = conf.get("tool_pair", "metaphlan,bracken")
    out_dir = conf.get("out_dir", os.path.join(conf.get("base_out_dir", "diversity_results"), "eco_metrics"))
    log_level = conf.get("log_level", "INFO")

    agreement_mode = conf.get("agreement_mode", "clr")
    subcomp_mode = conf.get("subcomp_mode", "prevalence")
    subcomp_prevalence = float(conf.get("subcomp_prevalence", 0.5))
    subcomp_topk = int(conf.get("subcomp_topk", 0))
    zero_replacement = conf.get("zero_replacement", "pseudocount")
    pseudocount = float(conf.get("pseudocount", 1e-6))
    zero_delta = float(conf.get("zero_delta", 1e-6))
    dirichlet_alpha = float(conf.get("dirichlet_alpha", 0.5))
    bias_correct = conf.get("bias_correct", "none")
    bias_shrink = float(conf.get("bias_shrink", 0.1))
    deming_delta = float(conf.get("deming_delta", 1.0))
    detect_thresh = conf.get("detect_thresh", "fixed")
    detect_alpha = float(conf.get("detect_alpha", 0.01))
    detect_fixed_t = int(conf.get("detect_fixed_t", 1))
    target_detect = float(conf.get("target_detect", 0.95))
    depths = conf.get("depths")
    forecast_depths = conf.get("forecast_depths")
    beta_grid = conf.get("beta_grid")
    metadata = conf.get("metadata")
    sample_col = conf.get("sample_col")
    include_regex = conf.get("include_regex")
    cohort_col = conf.get("cohort_col")
    cohort_allow = conf.get("cohort_allow")
    date_col = conf.get("date_col")
    exclude_old_before = conf.get("exclude_old_before")
    neg_regex = conf.get("neg_regex", "(?i)(neg|blank|ntc|control)")
    n_bootstrap = int(conf.get("n_bootstrap", 0))
    min_occupancy = float(conf.get("min_occupancy", 0.5))
    seed = int(conf.get("seed", 1337))

    if not merged_path or not counts_path or not rank or not out_dir:
        raise ValueError("Missing required ecology inputs: merged, counts, rank, out_dir")

    logger = _eco_setup_logging(log_level)
    _eco_ensure_dir(out_dir)

    merged = pd.read_csv(merged_path, sep="\t")
    t1, t2 = [s.strip() for s in tool_pair.split(",")] if "," in tool_pair else ("metaphlan","bracken")
    for col in ["sample","taxon","rank",t1,t2]:
        if col not in merged.columns:
            raise ValueError(f"Missing column in merged: {col}")
    if "rank" in merged.columns:
        merged = merged[merged["rank"] == rank].copy()
    merged = merged.sort_values(["sample","taxon"], kind="mergesort").reset_index(drop=True)

    pos_samples, neg_samples, _exc = _eco_build_sample_lists(
        merged=merged,
        metadata_path=metadata,
        sample_col=sample_col,
        include_regex=include_regex,
        cohort_col=cohort_col,
        cohort_allow=cohort_allow,
        date_col=date_col,
        exclude_old_before=exclude_old_before,
        neg_regex=neg_regex,
        out_dir=out_dir,
        logger=logger,
    )
    pos_df = merged[merged["sample"].isin(pos_samples)].copy()

    keep_taxa = _eco_build_fixed_subcomposition(pos_df, pos_samples, t1, t2,
                                                mode=subcomp_mode,
                                                prevalence=subcomp_prevalence,
                                                topk=subcomp_topk)
    pd.DataFrame({"taxon": keep_taxa}).to_csv(os.path.join(out_dir, "subcomposition_taxa.tsv"), sep="\t", index=False)

    bias_df = None
    if bias_correct == "clr-mean" and agreement_mode == "clr":
        bias_df = _eco_estimate_clr_bias(merged, pos_samples, keep_taxa, t1, t2,
                                         zero_replacement, pseudocount, zero_delta, dirichlet_alpha,
                                         lam=bias_shrink)
        bias_df.to_csv(os.path.join(out_dir, "tool_bias.tsv"), sep="\t", index=False)

    rng_main = np.random.default_rng(int(seed))
    agree = _eco_agreement_per_sample_fixed(
        pos_df, t1, t2, keep_taxa,
        mode=agreement_mode,
        zero_replacement=zero_replacement, pseudocount=pseudocount,
        alr_ref=conf.get("alr_ref"),
        alr_topk=int(conf.get("alr_topk", 0)),
        plr_maxpairs=int(conf.get("plr_maxpairs", 0)),
        rng=rng_main,
        zero_delta=zero_delta, dirichlet_alpha=dirichlet_alpha,
        bias_df=bias_df, bias_mode=bias_correct, deming_delta=deming_delta
    )
    agree.to_csv(os.path.join(out_dir, f"agreement_per_sample_{agreement_mode}.tsv"), sep="\t", index=False)

    summary: List[str] = []
    def _add_sum(line: str): summary.append(line)
    def _safe_med(x):
        a = pd.to_numeric(x, errors='coerce').astype(float).values
        m = np.isfinite(a)
        return float(np.nanmedian(a[m])) if m.any() else float('nan')
    _add_sum(f"samples_pos={len(set(pos_samples))} samples_neg={len(set(neg_samples))} taxa={pos_df['taxon'].nunique()} rank={rank}")
    _add_sum(f"Agreement ({agreement_mode}) median: pearson={_safe_med(agree['pearson']):.3f} spearman={_safe_med(agree['spearman']):.3f} aitchison={_safe_med(agree['aitchison_clr_norm']):.3f}")

    counts_df = None
    if counts_path and os.path.exists(counts_path):
        counts_df = pd.read_csv(counts_path, sep="\t" if counts_path.lower().endswith((".tsv",".txt")) else ",")[
            ["sample","N"]]
        counts_df = counts_df[counts_df["sample"].isin(pos_samples + neg_samples)].copy()
    else:
        logger.warning("No counts file provided; skipping moments/AFD/occupancy/richness/forecasts/GOF/rarefaction.")

    if counts_df is not None and not counts_df.empty:
        counts = _eco_build_counts_from_pct(merged[merged["sample"].isin(pos_samples + neg_samples)][["sample","taxon",t2]].copy(), counts_df, t2)

        # thresholds (kFDR/neg-poisson) omitted in trimmed version; use simple n>0
        counts_pos = counts[counts["sample"].isin(pos_samples)].copy()
        Ns_pos = counts_df[counts_df["sample"].isin(pos_samples)].drop_duplicates("sample")["N"].astype(float).values
        medN = float(np.median(Ns_pos)) if Ns_pos.size>0 else np.nan

        mom = _eco_unbiased_moments(counts_pos)

        bar_x = mom.set_index("taxon")["mu1"].astype(float)
        beta_grid_vals = _eco_parse_csv_list(beta_grid)
        occ_emp = counts_pos.groupby("taxon")["n"].apply(lambda x: float((x>0).mean()))
        beta_shared = _eco_fit_shared_beta(bar_x.values, Ns_pos, occ_emp.reindex(bar_x.index).fillna(0.0).values, beta_grid_vals)

        gof = _eco_chisq_gof_nb(counts_pos, bar_x, beta_shared)
        if gof is None or gof.empty:
            gof = pd.DataFrame({"taxon": [], "gof_chisq": [], "gof_df": [], "gof_p": []})
        else:
            gof["gof_q"] = _eco_bh_fdr(gof["gof_p"].values)
            gof["gof_pass"] = (gof["gof_q"] <= 0.05)

        rows = []
        for tx in bar_x.index:
            bx = float(bar_x.loc[tx]); b = float(beta_shared)
            occ_cur = float(occ_emp.get(tx, np.nan))
            Nstar = _eco_n_star_q(bx, b if np.isfinite(b) else np.inf, target_detect) if bx>0 else np.nan
            rows.append({"taxon": tx, "bar_x": bx, "beta": float(b), "occ_current_mean": occ_cur, "N_star_q": Nstar, "theta_hat": np.nan})
        afd = pd.DataFrame(rows).sort_values("taxon", kind="mergesort").reset_index(drop=True)
        if gof is not None and not gof.empty:
            afd = afd.merge(gof, on="taxon", how="left")
        afd.to_csv(os.path.join(out_dir, "afd_fit.tsv"), sep="\t", index=False)

        pred_occ = _eco_occupancy_from_beta(bar_x.values, Ns_pos, float(np.nanmedian(afd["beta"].astype(float))))
        emp_occ_vals = afd.set_index("taxon")["occ_current_mean"].reindex(bar_x.index).fillna(0.0).values
        rmse = float(np.sqrt(np.mean((pred_occ - emp_occ_vals)**2)))
        r = _eco_safe_corr(pred_occ, emp_occ_vals)
        r2 = (r*r) if np.isfinite(r) else np.nan
        _add_sum(f"Occupancy calibration: RMSE={rmse:.3f} R2={r2:.3f}")

        depth_rows = []
        for _, r in afd.iterrows():
            b = float(r["beta"]) if np.isfinite(r["beta"]) else np.inf
            bx = float(r["bar_x"])
            Nstar = float(r["N_star_q"]) if np.isfinite(r["N_star_q"]) else np.nan
            occ_medN = 1.0 - _eco_safe_nb_zero_prob(bx * medN, b) if np.isfinite(medN) else np.nan
            depth_rows.append({
                "taxon": r["taxon"], "bar_x": bx, "beta": b, "N_star_q": Nstar,
                "medianN": medN, "is_sufficient_at_current": bool(np.isfinite(Nstar) and medN >= Nstar),
                "occupancy_at_medianN": occ_medN, "observed_occupancy": float(occ_emp.get(r["taxon"], np.nan)),
                "theta_hat": r.get("theta_hat", np.nan)
            })
        depth_suff = pd.DataFrame(depth_rows).sort_values("taxon", kind="mergesort")
        depth_suff.to_csv(os.path.join(out_dir, "depth_sufficiency_per_species.tsv"), sep="\t", index=False)
        suff_pct = float(np.nanmean(depth_suff["is_sufficient_at_current"])) * 100.0
        _add_sum(f"Depth sufficiency at medianN: {suff_pct:.1f}% species; medianN={medN:.3g}")

        rich_rows = []
        dvals = _eco_parse_csv_list(depths) or []
        if medN > 0 and (not dvals or all(abs(d - medN) > 1e-9 for d in dvals)):
            dvals = [medN] + dvals
        beta_for_R = float(np.nanmedian(afd["beta"].astype(float)))
        bx_vec = bar_x.values
        for d in dvals:
            R, Rp = _eco_richness_and_slope(bx_vec, beta_for_R, float(d))
            rich_rows.append({"stratum": "ALL", "depth": float(d), "R": R, "R_prime": Rp})
        rich_df = pd.DataFrame(rich_rows)
        if not rich_df.empty:
            rich_df.sort_values(["stratum","depth"], kind="mergesort").to_csv(os.path.join(out_dir, "richness_vs_depth.tsv"), sep="\t", index=False)

        fdepths = _eco_parse_csv_list(forecast_depths) or []
        if fdepths:
            beta_used = beta_for_R
            rows2 = []
            for d in fdepths:
                for s, g in pos_df.groupby("sample", sort=False):
                    spp = g["taxon"].unique()
                    bx = bar_x.reindex(spp).dropna().values
                    if bx.size == 0:
                        zf = 1.0
                    else:
                        zf = float(np.mean(_eco_safe_nb_zero_prob(bx * float(d), beta_used)))
                    sugg = "pearson" if zf < 0.2 else "spearman"
                    rows2.append({"depth": float(d), "sample": s, "pred_zero_fraction": zf, "suggested": sugg})
            if rows2:
                fdf = pd.DataFrame(rows2).sort_values(["depth","sample"], kind="mergesort")
                fdf.to_csv(os.path.join(out_dir, "agreement_forecast_vs_depth.tsv"), sep="\t", index=False)
                agg = fdf.groupby("depth")["suggested"].value_counts(normalize=True).rename("prop").reset_index()
                agg_pivot = agg.pivot(index="depth", columns="suggested", values="prop").fillna(0.0).reset_index()
                agg_pivot.to_csv(os.path.join(out_dir, "agreement_forecast_summary.tsv"), sep="\t", index=False)
                for _, r in agg_pivot.iterrows():
                    _add_sum(f"Forecast at depth={r['depth']:.3g}: pearson%={100*r.get('pearson',0):.1f} spearman%={100*r.get('spearman',0):.1f}")

        rgrid = _eco_parse_csv_list(conf.get("rarefy_grid"))
        if rgrid:
            rare = _eco_rarefy_counts(counts_pos, [int(d) for d in rgrid], draws=50, rng=rng_main)
            if not rare.empty:
                agg = rare.groupby(["sample","depth"])["richness"].agg(["mean","median","std"]).reset_index()
                agg.to_csv(os.path.join(out_dir, "rarefaction_empirical.tsv"), sep="\t", index=False)

        pc_grid = _eco_parse_csv_list(conf.get("pseudocount_grid"))
        if pc_grid:
            rows3=[]
            for pc in pc_grid:
                a = _eco_agreement_per_sample_fixed(
                    pos_df, t1, t2, keep_taxa,
                    mode=agreement_mode,
                    zero_replacement="pseudocount", pseudocount=float(pc),
                    alr_ref=conf.get("alr_ref"),
                    alr_topk=int(conf.get("alr_topk", 0)),
                    plr_maxpairs=int(conf.get("plr_maxpairs", 0)),
                    rng=rng_main,
                    zero_delta=zero_delta, dirichlet_alpha=dirichlet_alpha,
                    bias_df=bias_df, bias_mode=bias_correct, deming_delta=deming_delta
                )
                rows3.append({"pseudocount": pc, "pearson_med": float(np.nanmedian(a["pearson"])), "spearman_med": float(np.nanmedian(a["spearman"]))})
            pd.DataFrame(rows3).to_csv(os.path.join(out_dir, f"agreement_sensitivity_{agreement_mode}.tsv"), sep="\t", index=False)

        qg = _eco_parse_csv_list(conf.get("q_grid"))
        if qg:
            rows4=[]
            for q in qg:
                for _, r in afd.iterrows():
                    bx = float(r["bar_x"]); b = float(r["beta"]) if np.isfinite(r["beta"]) else np.inf
                    Nstar = _eco_n_star_q(bx, b, q) if (np.isfinite(bx) and bx>0) else np.nan
                    rows4.append({"taxon": r["taxon"], "q": q, "N_star_q": Nstar})
            pd.DataFrame(rows4).to_csv(os.path.join(out_dir, "nstar_grid.tsv"), sep="\t", index=False)

    with open(os.path.join(out_dir, "summary.txt"), "w") as f:
        for line in summary:
            f.write(line + "\n")
    return 0


if __name__ == "__main__":
    main()
