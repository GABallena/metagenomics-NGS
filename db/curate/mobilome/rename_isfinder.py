# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO
import csv

# Define input and output
input_fasta = Path("ISFinder_db/ISfinder_sequences.fasta")
output_fasta = Path("ISFinder_db_renamed.fasta")
output_metadata = Path("ISFinder_db_metadata.tsv")

renamed_records = []
metadata_rows = []

# Counter
element_counter = 1

# Process each record
for record in SeqIO.parse(input_fasta, "fasta"):
    original_id = record.id
    new_id = f"ISFinder_IS_Element{element_counter}"
    element_counter += 1

    record.id = new_id
    record.description = ""
    renamed_records.append(record)

    metadata_rows.append([
        new_id,
        "ISFinder",
        "Insertion_Sequence",
        original_id
    ])

# Write outputs
SeqIO.write(renamed_records, output_fasta, "fasta")

with open(output_metadata, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ShortBRED_ID", "Tool", "Category", "Original_ID"])
    writer.writerows(metadata_rows)
