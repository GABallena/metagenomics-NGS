## Experimental Workflows

Exploratory scripts and Snakemake DAGs used during pipeline prototyping. They are kept for transparency and reproducibility, but they are not part of the production workflow.

- `trim_randomizer.smk`: parameter sweep generator for multiple trimmers. Produces synthetic TSV logs of the random seeds/params used in each iteration.
- `trimming_qc_optimization.py`: standalone Trimmomatic + FastQC harness for quick local trials; expects `raw_reads/` populated by the user.
- `raw_vs_trimmed_analysis.py`: exploratory statistics on FASTQ pairs (length, GC, entropy, clustering). Safe to run without data present; it emits a friendly message if inputs are missing.

Guidance:
- Run inside an isolated environment; no data is tracked in the repo.
- Treat outputs as disposable sandboxes; prefer `workflows/` for reusable modules.
