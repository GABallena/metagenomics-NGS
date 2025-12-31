#!/bin/bash
# ==========================
# Batch download ARG DBs
# ==========================
# Requires: aria2c, unzip/tar

set -euo pipefail

DB_ROOT="${DB_ROOT:-db/downloads/amr}"
mkdir -p "$DB_ROOT"

# Create a directory
mkdir -p ARG_DBs && cd ARG_DBs

# --------------------------
# CARD
# --------------------------
# Nucleotide FASTA (Protein Homolog Model Sequences)
aria2c -x 8 -s 8 -o card_nucleotide.fasta \
  "https://card.mcmaster.ca/latest/data/nucleotide_fasta_protein_homolog_model.fasta"

# Protein FASTA (optional)
aria2c -x 8 -s 8 -o card_protein.fasta \
  "https://card.mcmaster.ca/latest/data/protein_fasta_protein_homolog_model.fasta"

# Ontology (ARO term mappings)
aria2c -x 8 -s 8 -o card_ontology.json \
  "https://card.mcmaster.ca/latest/data/card.json"

# --------------------------
# AMRFinderPlus
# --------------------------
# NCBI FTP – nucleotide DB (latest)
aria2c -x 8 -s 8 -o amrfinder_nucl.fa.gz \
  "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/latest/amrfinder_nucl.fa.gz"

# Protein DB
aria2c -x 8 -s 8 -o amrfinder_prot.fa.gz \
  "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/latest/amrfinder_prot.fa.gz"

# --------------------------
# SARG (Structured ARG DB)
# --------------------------
aria2c -x 8 -s 8 -o SARG.fasta \
  "http://smile.hku.hk/SARG/download/SARG.fasta"

# --------------------------
# ResFinder DB
# --------------------------
# Latest ResFinder nucleotide database (master.zip from Bitbucket)
aria2c -x 8 -s 8 -o resfinder_db.zip \
  "https://bitbucket.org/genomicepidemiology/resfinder_db/get/master.zip"
unzip resfinder_db.zip -d resfinder_db

# --------------------------
# MEGARes v3.0
# --------------------------
aria2c -x 8 -s 8 -o megares_v3.0.fa \
  "https://megares.meglab.org/download/megares_database_v3.0.fasta"

# Annotation file (class/ontology)
aria2c -x 8 -s 8 -o megares_annotations_v3.0.csv \
  "https://megares.meglab.org/download/megares_annotations_v3.0.csv"

# --------------------------
# DeepARG DB
# --------------------------
# From Zenodo release
aria2c -x 8 -s 8 -o deeparg_db.zip \
  "https://zenodo.org/records/8280582/files/deeparg.zip"
unzip deeparg_db.zip -d deeparg_db

# ==========================
echo "All downloads attempted. Verify file integrity (md5/sha256)!"
