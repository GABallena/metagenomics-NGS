# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO
import csv

# Define the source directory and output paths
input_dir = Path("Integrall_db")
output_fasta = Path("Integrall_db_renamed.fasta")
output_metadata = Path("Integrall_db_metadata.tsv")

# Prepare output files
renamed_records = []
metadata_rows = []

# Counter to assign unique IDs
gene_counter = 1

# Iterate over all FASTA files in the directory
for fasta_file in sorted(input_dir.glob("*.fasta")):
    for record in SeqIO.parse(fasta_file, "fasta"):
        original_id = record.id
        filename_stem = fasta_file.stem

        # Create new ID
        new_id = f"Integrall_Integron_Gene{gene_counter}"
        gene_counter += 1

        # Update FASTA record
        record.id = new_id
        record.description = ""
        renamed_records.append(record)

        # Append metadata
        metadata_rows.append([
            new_id,
            "Integrall",
            "Integron",
            original_id,
            filename_stem
        ])

# Write combined and renamed FASTA
SeqIO.write(renamed_records, output_fasta, "fasta")

# Write metadata table
with open(output_metadata, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["ShortBRED_ID", "Tool", "Category", "Original_ID", "Source_File"])
    writer.writerows(metadata_rows)
