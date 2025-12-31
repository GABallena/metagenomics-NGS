# BUSCO SCG Workflow (Anonymized)

This note documents a cleaned, shareable version of the BUSCO Single Copy Gene (SCG) workflow used to seed k-mer databases and accelerate BLAST/core database lookups. All identifiers and paths are generic so you can adapt them to any environment.

## Objectives
- Enumerate SCG datasets from BUSCO metadata.
- Optionally download BUSCO lineages and the `core_nt` database.
- Extract ancestral/variant FASTA files and generate k-mer databases.
- Provide lightweight utilities to summarize BLAST hits.

## Workflow Outline
1. **List BUSCO datasets**  
   ```bash
   busco --list-datasets
   busco download --all   # optional bulk download
   ```

2. **Extract SCG dataset IDs**  
   ```python
   # file_versions_busco.py
   import pandas as pd

   df = pd.read_csv("busco_downloads/file_versions.tsv", sep="\t", header=None)
   df.columns = ["Dataset", "Date", "Hash", "Domain", "Type"]
   df["Dataset"].to_csv("all_scg_list.txt", index=False, header=False)
   ```

3. **(Optional) Download `core_nt` with aria2c**  
   ```bash
   # download_core_nt.sh
   for i in {00..68}; do
     aria2c -x 16 "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/core_nt.${i}.tar.gz"
   done
   ```

4. **Extract ancestral FASTAs from BUSCO lineages**  
   ```python
   # extract_ancestral_and_variants.py
   import os, shutil, pathlib

   src = pathlib.Path("busco_downloads/lineages")
   dest = pathlib.Path("extracted_fastas"); dest.mkdir(exist_ok=True)

   for root, _, files in os.walk(src):
       for name in files:
           if name in {"ancestral", "ancestral_variants"}:
               lineage = pathlib.Path(root).name
               out = dest / f"{lineage}_{name}.fasta"
               shutil.copy(pathlib.Path(root) / name, out)
   ```

5. **Generate k-mer databases (k=21..100)**  
   ```python
   # generate_kmers.py
   import os

   def kmers(seq, k): return {seq[i:i+k] for i in range(len(seq)-k+1)}

   os.makedirs("kmer_databases", exist_ok=True)
   for fasta in os.listdir("extracted_fastas"):
       if not fasta.endswith(".fasta"): continue
       seq = "".join(
           line.strip() for line in open(os.path.join("extracted_fastas", fasta))
           if not line.startswith(">")
       )
       lineage = fasta.split("_")[0]
       for k in range(21, 101):
           out = f"kmer_databases/{lineage}_k{k}.txt"
           with open(out, "w") as handle:
               handle.write("\n".join(kmers(seq, k)))
   ```

6. **Summarize BLAST hits**  
   ```python
   # nr_busco.py
   import pandas as pd

   columns = ["QueryID", "SubjectID", "Identity", "Length", "Evalue"]
   df = pd.read_csv("blast_results.txt", sep="\t", header=None, names=columns)
   df.sort_values("Evalue").drop_duplicates("QueryID").to_csv("top_hits.csv", index=False)
   ```

## Tips for Adapting
- Keep downloads and generated artifacts under `data/` (not tracked).
- Replace hardcoded paths with environment variables when deploying on HPC.
- Consider Snakemake wrappers for each step if you want this in a larger DAG.

## Provenance
This document is anonymized and designed for portfolio use; no proprietary data or identifiers are included.
