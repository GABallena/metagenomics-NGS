import pandas as pd

# Load BLAST results
results_file = "scg_vs_nr_results.txt"
columns = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
blast_results = pd.read_csv(results_file, sep="\t", names=columns)

# Filter top hits per query
top_hits = blast_results.sort_values(["qseqid", "evalue"]).drop_duplicates("qseqid")

# Save filtered results
top_hits.to_csv("top_scg_hits.csv", index=False)
print("Filtered results saved to 'top_scg_hits.csv'.")

