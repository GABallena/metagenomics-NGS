# Portfolio-safe copy: paths/identifiers generalized; databases not included.
from pathlib import Path
from Bio import SeqIO, Seq
import csv
import re

# === Input/Output Paths ===
input_fasta = Path("mobileOGdb_90-2.0.fasta")
output_fasta = Path("mobileOGdb_renamed.fasta")
output_metadata = Path("mobileOGdb_metadata.tsv")

# === Setup ===
renamed_records = []
metadata_rows = []

def clean_sequence(seq):
    """Removes invalid characters from a sequence, keeping only ACGTN."""
    return re.sub(r"[^ACGTNacgtn]", "", str(seq))

# === Process FASTA Records ===
for i, record in enumerate(SeqIO.parse(input_fasta, "fasta"), 1):
    original_id = record.id
    new_id = f"MobileOG_Entry{i}"
    
    # Clean and update
    cleaned_seq = clean_sequence(record.seq)
    record.seq = Seq.Seq(cleaned_seq)
    record.id = new_id
    record.description = ""

    renamed_records.append(record)
    metadata_rows.append([new_id, "MobileOG_DB", "MGE_Annotation", original_id])

# === Write outputs ===
SeqIO.write(renamed_records, output_fasta, "fasta")

with open(output_metadata, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["MobileOG_ID", "Tool", "Category", "Original_ID"])
    writer.writerows(metadata_rows)

print(f"✅ Done: {len(renamed_records)} records written to {output_fasta.name}")
