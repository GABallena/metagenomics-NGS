# Database utilities (fetch → curate → harmonize)

This `db/` module contains **scripts only** for acquiring and preparing reference databases used in AMR / mobilome / pathogen-focused metagenomics workflows.

## Golden rules (public-safe)
- Do **not** commit downloaded FASTA/GenBank/TSV databases, genome archives, or large outputs.
- Each upstream database has its own terms. Follow licensing/ToS for download and use.
- Commit scripts + documentation only. Keep outputs outside the repo or in ignored folders.

Suggested `.gitignore` additions (repo root):

```gitignore
# downloaded databases / genomes / outputs
db/**/downloads/
db/**/work/
db/**/*.fasta
db/**/*.fa
db/**/*.fna
db/**/*.ffn
db/**/*.fas
db/**/*.fsa
db/**/*.gb
db/**/*.gbk
db/**/*.gz
db/**/*.zip
db/**/*_metadata.tsv
db/**/*_renamed*
reference_genomes/
```

## Layout
- `db/fetch/` acquisition helpers (NCBI pulls, genome fetchers)
- `db/curate/` source-specific renamers + metadata emitters
- `db/harmonize/` cross-database cleanup/mapping utilities

## Typical flow
1) Fetch raw sources locally  
2) Curate into stable IDs + metadata TSV per source  
3) Harmonize / map across sources as needed  

See subfolder READMEs for notes and prerequisites.
