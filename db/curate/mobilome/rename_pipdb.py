# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO
import csv

# Input/output paths
input_fasta = Path("PIP_db/PIP_db.fasta")
output_fasta = Path("PIP_db_renamed.fasta")
output_metadata = Path("PIP_db_metadata.tsv")

renamed_records = []
metadata_rows = []

plasmid_counter = 1

for record in SeqIO.parse(input_fasta, "fasta"):
    original_id = record.id
    new_id = f"PIPdb_Pathogenic_Plasmid{plasmid_counter}"
    plasmid_counter += 1

    record.id = new_id
    record.description = ""
    renamed_records.append(record)

    metadata_rows.append([
        new_id,
        "PIPdb",
        "Pathogenic_Plasmid",
        original_id
    ])

# Write renamed FASTA
SeqIO.write(renamed_records, output_fasta, "fasta")

# Write metadata TSV
with open(output_metadata, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ShortBRED_ID", "Tool", "Category", "Original_ID"])
    writer.writerows(metadata_rows)
