# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO
import csv

# Input/output
input_dir = Path("plasmidfinder_db")
output_fasta = Path("plasmidfinder_db_renamed.fasta")
output_metadata = Path("plasmidfinder_db_metadata.tsv")

renamed_records = []
metadata_rows = []
replicon_counter = 1

# Loop through all .fsa files
for fasta_file in sorted(input_dir.glob("*.fsa")):
    for record in SeqIO.parse(fasta_file, "fasta"):
        original_id = record.id
        new_id = f"PlasmidFinder_Replicon{replicon_counter}"
        replicon_counter += 1

        record.id = new_id
        record.description = ""
        renamed_records.append(record)

        metadata_rows.append([
            new_id,
            "PlasmidFinder",
            "Replicon",
            original_id,
            fasta_file.name
        ])

# Write outputs
SeqIO.write(renamed_records, output_fasta, "fasta")

with open(output_metadata, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ShortBRED_ID", "Tool", "Category", "Original_ID", "Source_File"])
    writer.writerows(metadata_rows)
