# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO
import csv

# Define paths
input_fasta = Path("IMG_PR_db/IMGPR_nucl.fna")
output_fasta = Path("IMGPR_db_renamed.fasta")
output_metadata = Path("IMGPR_db_metadata.tsv")

# Prepare containers
renamed_records = []
metadata_rows = []

# Counter for unique IDs
gene_counter = 1

# Parse and rename
for record in SeqIO.parse(input_fasta, "fasta"):
    original_id = record.id
    new_id = f"IMGPR_Plasmid_Gene{gene_counter}"
    gene_counter += 1

    # Update record
    record.id = new_id
    record.description = ""
    renamed_records.append(record)

    # Add metadata row
    metadata_rows.append([
        new_id,
        "IMGPR",
        "Plasmid",
        original_id
    ])

# Write renamed FASTA
SeqIO.write(renamed_records, output_fasta, "fasta")

# Write metadata
with open(output_metadata, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ShortBRED_ID", "Tool", "Category", "Original_ID"])
    writer.writerows(metadata_rows)
