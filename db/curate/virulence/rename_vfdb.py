# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO
import csv

# Paths
input_fasta = Path("VF_db/VFDB_setA_nt.fas")
output_fasta = Path("VF_db_renamed.fasta")
output_metadata = Path("VF_db_metadata.tsv")

renamed_records = []
metadata_rows = []

vf_counter = 1

for record in SeqIO.parse(input_fasta, "fasta"):
    original_id = record.id
    new_id = f"VFDB_VirulenceFactor{vf_counter}"
    vf_counter += 1

    record.id = new_id
    record.description = ""
    renamed_records.append(record)

    metadata_rows.append([
        new_id,
        "VFDB",
        "Virulence_Factor",
        original_id
    ])

# Write renamed FASTA
SeqIO.write(renamed_records, output_fasta, "fasta")

# Write metadata
with open(output_metadata, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ShortBRED_ID", "Tool", "Category", "Original_ID"])
    writer.writerows(metadata_rows)
