# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO
import csv

# Paths
input_fasta = Path("TnCentral_db/tncentral_integrall_isfinder.fa")
output_fasta = Path("TnCentral_db_renamed.fasta")
output_metadata = Path("TnCentral_db_metadata.tsv")

renamed_records = []
metadata_rows = []

mge_counter = 1

for record in SeqIO.parse(input_fasta, "fasta"):
    original_id = record.id
    new_id = f"TnCentral_MGE{mge_counter}"
    mge_counter += 1

    record.id = new_id
    record.description = ""
    renamed_records.append(record)

    metadata_rows.append([
        new_id,
        "TnCentral",
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
