# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO
import csv

# Define input/output paths
input_fasta = Path("MGEfinder_db/mgedb/data/sequences.d/mge_records.fna")
output_fasta = Path("MGEfinder_db_renamed.fasta")
output_metadata = Path("MGEfinder_db_metadata.tsv")

renamed_records = []
metadata_rows = []

element_counter = 1

# Parse and rename
for record in SeqIO.parse(input_fasta, "fasta"):
    original_id = record.id
    new_id = f"MGEfinder_mobile_element{element_counter}"
    element_counter += 1

    # Update record
    record.id = new_id
    record.description = ""
    renamed_records.append(record)

    # Metadata row
    metadata_rows.append([
        new_id,
        "MGEfinder",
        "Mobile_Genetic_Element",
        original_id
    ])

# Write renamed FASTA
SeqIO.write(renamed_records, output_fasta, "fasta")

# Write metadata TSV
with open(output_metadata, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ShortBRED_ID", "Tool", "Category", "Original_ID"])
    writer.writerows(metadata_rows)
