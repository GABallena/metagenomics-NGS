# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO
import csv

# Paths
input_fasta = Path("TREP_db/TREP_db.fasta")
output_fasta = Path("TREP_db_renamed.fasta")
output_metadata = Path("TREP_db_metadata.tsv")

renamed_records = []
metadata_rows = []

te_counter = 1

for record in SeqIO.parse(input_fasta, "fasta"):
    original_id = record.id
    new_id = f"TREP_TransposableElement{te_counter}"
    te_counter += 1

    record.id = new_id
    record.description = ""
    renamed_records.append(record)

    metadata_rows.append([
        new_id,
        "TREP",
        "Transposable_Element",
        original_id
    ])

# Write renamed FASTA
SeqIO.write(renamed_records, output_fasta, "fasta")

# Write metadata TSV
with open(output_metadata, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ShortBRED_ID", "Tool", "Category", "Original_ID"])
    writer.writerows(metadata_rows)
