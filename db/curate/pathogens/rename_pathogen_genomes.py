# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO
import csv
from collections import defaultdict

# Directories
input_dir = Path("pathogen_db")
output_fasta = Path("pathogen_db_renamed.fasta")
output_metadata = Path("pathogen_db_metadata.tsv")

renamed_records = []
metadata_rows = []
species_counters = defaultdict(int)

# Go through species folders
for species_dir in sorted(input_dir.iterdir()):
    if not species_dir.is_dir():
        continue
    species_name = species_dir.name.replace(" ", "_")
    fna_files = list(species_dir.glob("**/GCF*/*_genomic.fna"))
    for fna in fna_files:
        species_counters[species_name] += 1
        genome_id = f"{species_name}_genome{species_counters[species_name]}"
        for record in SeqIO.parse(fna, "fasta"):
            original_id = record.id
            record.id = genome_id
            record.description = ""
            renamed_records.append(record)

            metadata_rows.append([
                genome_id,
                "NCBI_Datasets",
                "Pathogen_Genome",
                original_id,
                species_name,
                fna.name
            ])

# Write combined FASTA
SeqIO.write(renamed_records, output_fasta, "fasta")

# Write metadata
with open(output_metadata, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Pathogen_ID", "Tool", "Category", "Original_ID", "Species", "Source_File"])
    writer.writerows(metadata_rows)
