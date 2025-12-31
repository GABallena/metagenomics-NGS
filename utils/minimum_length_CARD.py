from pathlib import Path
import argparse
from Bio import SeqIO


def find_min_length(fasta_file: Path) -> int:
    min_length = float("inf")
    for record in SeqIO.parse(str(fasta_file), "fasta"):
        seq_length = len(record.seq)
        if seq_length < min_length:
            min_length = seq_length
    return int(min_length)


def main():
    parser = argparse.ArgumentParser(
        description="Report minimum sequence length in a FASTA file (e.g., CARD)."
    )
    parser.add_argument(
        "fasta",
        nargs="?",
        default=Path("test/databases/CARD_sequences/extracted/nucleotide_fasta_protein_homolog_model.fasta"),
        type=Path,
        help="Path to FASTA file (default: relative test database location).",
    )
    args = parser.parse_args()

    if not args.fasta.exists():
        raise FileNotFoundError(f"FASTA file not found: {args.fasta}")

    min_length = find_min_length(args.fasta)
    print(f"Minimum sequence length in {args.fasta}: {min_length} bp")


if __name__ == "__main__":
    main()
