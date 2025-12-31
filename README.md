# Metagenomics NGS Portfolio Workflow

Metagenomics processing toolkit with production-grade Snakemake workflows, HPC-ready batch scripts, and exploratory notebooks/scripts for QC and downstream analysis. All identifiers are anonymized and paths are relative so the project can be shared as a portfolio example or used as a clean starting point.

## What's Inside
- `metagenomics_general.smk`: end-to-end short-read workflow (trimming, FastQC, Kraken2, Bracken).
- `config/config.yaml`: single source of truth for inputs, directories, parameters, threads, and database locations.
- `profiles/slurm/`: portable Slurm templates for running the pipeline and staging reads.
- `analysis/`, `experiments/`, `workflows/`: curated R/Python/Snakemake prototypes for QC, parameter sweeps, and modular sub-workflows.
- `docs/`: anonymized documentation describing SCG/BUSCO workflows and dataset preparation.
- `tools/` and `utils/`: helper scripts for scaffolding and data transfer.

## Quickstart
1) Install Snakemake (via conda or mamba).  
2) Edit `config/config.yaml` with your read paths, adapter file, and Kraken2 database location.  
3) Dry run to check the DAG:
   ```bash
   snakemake --snakefile metagenomics_general.smk -n
   ```
4) Execute locally (uses per-rule conda envs):
   ```bash
   snakemake --snakefile metagenomics_general.smk --use-conda --cores 8 --rerun-incomplete --printshellcmds
   ```
5) Run on Slurm:
   ```bash
   sbatch profiles/slurm/metagenomics_general.sbatch /path/to/project
   ```

## Project Layout
- `config/`: YAML configs for the main pipeline and modular workflows.
- `workflows/`: reusable Snakemake components (preprocessing, phylogenetics, MGE exploration).
- `experiments/`: parameter sweeps and exploratory scripts (kept for reproducibility, not production).
- `analysis/`: plotting and QC utilities.
- `profiles/`: cluster execution templates.
- `docs/`: workflow notes and references (anonymized).

## Notes on Anonymization
- No sample identifiers, institutions, or user-specific paths remain; everything is placeholder or relative.
- No raw data is bundled. Scripts assume you provide your own FASTQ/FASTA inputs.
- Replace demo paths in `config/config.yaml` before running on real data.

## Reproducibility
- All workflows are conda-backed; set `--conda-prefix` if you want an isolated env cache.
- For deterministic resource use, adjust `threads` and `mem_mb` in the config files.

## Contributing / Extending
- Add new tools by extending `config/params` and injecting rule params.
- Keep experimental code under `experiments/` and document design decisions in `docs/`.
- Use the Slurm templates as a base for other schedulers by swapping resource directives.
