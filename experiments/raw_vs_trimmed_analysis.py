#!/usr/bin/env python3
"""
Experimental analysis script for comparing raw vs trimmed FASTQ reads.

Purpose:
- Validate trimming effects on read quality, length, GC content
- Assess downstream impacts on sequence complexity and gene prediction
- Explore clustering behavior before and after trimming

Notes:
- Intended for exploratory / validation use, not production pipelines
- Assumes generic directory layout (raw_reads/, trimmed_reads/)
- No real data is included in this repository
"""

import gzip
import os
from pathlib import Path
import numpy as np
from collections import defaultdict, Counter
import math
import logging
from datetime import datetime
from tqdm import tqdm
import sys
import pandas as pd

try:
    import matplotlib.pyplot as plt
    PLOTTING_ENABLED = True
except ImportError:
    PLOTTING_ENABLED = False

from Bio import pairwise2
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.cluster import DBSCAN, KMeans, SpectralClustering
import networkx as nx


def setup_logging():
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"analysis_{timestamp}.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return log_file


def parse_quality_scores(qual_str):
    return [ord(c) - 33 for c in qual_str]


def calculate_gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100 if seq else 0


def calculate_shannon_entropy(sequence):
    counts = Counter(sequence)
    length = len(sequence)
    return -sum(
        (count / length) * math.log2(count / length)
        for count in counts.values()
        if length > 0
    )


def analyze_fastq_gz(filepath):
    logging.info(f"Analyzing FASTQ: {filepath}")

    stats = {
        "read_count": 0,
        "length_dist": defaultdict(int),
        "avg_qual_scores": [],
        "gc_content": [],
        "base_composition": defaultdict(int),
    }

    with gzip.open(filepath, "rt") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break

            seq = f.readline().strip()
            f.readline()
            qual = f.readline().strip()

            stats["read_count"] += 1
            stats["length_dist"][len(seq)] += 1
            stats["avg_qual_scores"].append(np.mean(parse_quality_scores(qual)))
            stats["gc_content"].append(calculate_gc_content(seq))

            for base in seq:
                stats["base_composition"][base] += 1

    return stats


def compare_raw_vs_trimmed():
    setup_logging()
    logging.info("Starting raw vs trimmed read comparison")

    raw_dir = Path("raw_reads")
    trimmed_dir = Path("trimmed_reads")

    results = {}

    for stage, directory in [("raw", raw_dir), ("trimmed", trimmed_dir)]:
        for fastq in directory.glob("*.fastq.gz"):
            sample = fastq.stem.split(".")[0]
            results.setdefault(sample, {})
            results[sample][stage] = analyze_fastq_gz(fastq)

    logging.info("Comparison complete")
    return results


if __name__ == "__main__":
    try:
        compare_raw_vs_trimmed()
    except FileNotFoundError:
        print("Expected directories 'raw_reads/' and 'trimmed_reads/' were not found.")
