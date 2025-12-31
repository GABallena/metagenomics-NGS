# metagenomics-NGS

Modular metagenomics/NGS workflows and helper utilities for reproducible QC → assembly → binning → profiling → AMR/mobilome analyses, with HPC/Slurm-ready templates.

This repository is organized as:
- Snakemake workflows (`workflows/`, `metagenomics_general.smk`)
- Tool wrappers (`tools/`)
- Utilities (`utils/`)
- Analysis/plotting scripts (`analysis/`)
- Exploratory methods and one-off experiments (`experiments/`)
- Slurm execution profile + sbatch templates (`profiles/slurm/`)

## What this repo is (and is not)

### What it is
- A portfolio-grade collection of reproducible pipeline components.
- Focused on correct, automatable, scriptable analysis, designed to run on local Linux or HPC (Slurm).
- Modular: you can run a single workflow (e.g., preprocessing) or chain multiple stages.

### What it is not
- Not a “one-click turnkey pipeline” for every dataset.
- Does not ship any proprietary databases or sample data.
- Not tied to any institution or private project; paths and identifiers are generalized.

## Repository structure

```
.
├── metagenomics_general.smk           # orchestrator (high-level)
├── workflows/                         # canonical Snakemake modules
│   ├── preprocessing.smk
│   ├── assembly.smk
│   ├── assembly-dedup.smk
│   ├── binning.smk
│   ├── ARG.smk
│   ├── MGE.smk
│   ├── phylogenetics.smk
│   ├── config/
│   │   ├── preprocessing.yaml
│   │   └── phylogenetics.yaml
│   └── scripts/
│       ├── calculate_contig_quality.py
│       └── calculate_plasmid_percentage.py
├── profiles/slurm/                    # snakemake --profile slurm
│   ├── README.md
│   ├── download_k2db.sbatch
│   ├── metagenomics_general.sbatch
│   └── move_reads.sbatch
├── tools/                             # wrappers around common tools
├── utils/                             # file plumbing, parsing, small helpers
├── analysis/                          # plotting and reporting scripts
├── experiments/                       # exploratory scripts & notes
└── config/
    └── config.yaml                    # general configuration entrypoint
```

## Typical workflow stages (high level)

Depending on which modules you run, the pipeline commonly follows:

1. **QC / preprocessing**
   - raw read QC (FastQC)
   - trimming (fastp or trimmomatic)
   - summary stats extraction

2. **Assembly**
   - MEGAHIT assembly
   - optional contig deduplication/scaffolding

3. **Binning / MAG QC**
   - binning modules (tooling depends on env)
   - contig quality metrics / BUSCO utilities

4. **Profiling**
   - taxonomic profiling (Kraken2 + Bracken, MetaPhlAn)
   - visualization support (Pavian)

5. **AMR and mobilome / MGE**
   - ARG workflows (e.g., CARD/NCBI-based workflows, filtering helpers)
   - mobilome analysis helpers

## Requirements

This repo assumes a Linux-like environment. Many scripts call external tools.

Minimum:
- `python >= 3.9`
- `snakemake` (7+ recommended)
- `conda` or `mamba` (recommended for env management)
- Standard CLI tools: `bash`, `awk`, `sed`, `gzip`, `unzip`

Common external tools used by scripts/workflows (not exhaustive):
- FastQC, fastp, Trimmomatic
- MEGAHIT, SPAdes
- Kraken2, Bracken
- MetaPhlAn
- Bowtie2
- BUSCO
- BLAST+
- seqtk
- jellyfish

Conda environment YAMLs live under `envs/` (intended as references/templates).

## Quickstart (local)

### 1) Install Snakemake
Example (conda):
```bash
conda create -n snakemake -c conda-forge -y snakemake
conda activate snakemake
```

### 2) Configure
Start by inspecting:
- `config/config.yaml`
- `workflows/config/preprocessing.yaml`
- `workflows/config/phylogenetics.yaml`

Edit paths for:
- input FASTQ locations
- output directories
- database paths (Kraken2, Bracken, etc.)
- tool parameters (threads, memory, read length, etc.)

### 3) Dry run
```bash
snakemake -n -p --snakefile metagenomics_general.smk
```

### 4) Run (local)
```bash
snakemake --snakefile metagenomics_general.smk --cores 8
```

## Running on HPC (Slurm)

A Slurm profile is provided under `profiles/slurm/`.

Typical usage:
```bash
snakemake \
  --snakefile metagenomics_general.smk \
  --profile profiles/slurm \
  --jobs 50
```

Also see `profiles/slurm/*.sbatch` for job submission templates and patterns.

## Notes on databases

This repository does not include databases.
You must download/build them locally and point configs/scripts to the correct paths.

Examples:
- Kraken2 DB: `tools/kraken_pipeline_script.sh` + `profiles/slurm/download_k2db.sbatch`
- Core nt downloads: `tools/downloadcore_nt.sh` (template for segmented download workflows)

## Reproducibility and safety

- Scripts are written to avoid hardcoded user paths and identity-bearing strings.
- Output directories are generally configurable via CLI args/env vars or Snakemake config.
- Do not commit raw data, intermediate BAMs, or large results into the repo.
- Prefer running with a clean directory structure and keep outputs outside the git working tree where possible.

## How to cite / reuse

If you reuse the workflow structure or scripts, a simple attribution in docs is appreciated.
No warranty: these are research/engineering utilities intended to be adapted to your environment.

## Contributing / Issues

If you find a bug:
- open an issue with minimal reproduction steps
- include tool versions and the relevant command used
- avoid sharing private sample identifiers or paths

## License
MIT
