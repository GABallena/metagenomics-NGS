import sys
import os
import zipfile
import tempfile
import csv

def parse_fastqc_summary(summary_file):
    """Parse the FastQC summary file and return a dictionary with status for each metric."""
    status = {}
    with open(summary_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 3:
                metric_key = shorten_metric(parts[1])
                status[metric_key] = parts[0]
    return status

def shorten_metric(metric_name):
    """Shorten the metric names to the provided abbreviations."""
    return {
        'Adapter Content': 'AC',
        'Per sequence base quality': 'BQ',
        'Per Base GC content': 'BGC',
        'Per sequence GC content': 'SGC',
        'Per base N content': 'SN',
        'Sequence Length Distribution': 'SL',
        'Sequence Duplication Levels': 'SD',
        'Basic Statistics': 'Stat',
        'Kmer Content': 'Kmer',
        'Per base sequence content': 'BSC',
        # Add other metrics if necessary
    }.get(metric_name, "UNK")  # Use "UNK" for unknown metrics

def extract_fastqc_zip(fastqc_zip):
    """Extract the FastQC zip file to a temporary directory and return the path to the extracted folder."""
    temp_dir = tempfile.mkdtemp()
    with zipfile.ZipFile(fastqc_zip, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)
    return temp_dir

def find_summary_file(extracted_dir):
    """Search for the summary.txt file in the extracted directory."""
    for root, dirs, files in os.walk(extracted_dir):
        if "summary.txt" in files:
            return os.path.join(root, "summary.txt")
    return None

def generate_report(fastqc_zips, output_tsv, legend_tsv):
    """Generates a TSV report with FastQC metrics for each sample."""
    all_metrics = set()

    # First pass: collect all unique metrics
    sample_statuses = {}
    for fastqc_zip in fastqc_zips:
        # Determine sample and read direction (R1 or R2)
        base_name = os.path.basename(fastqc_zip)
        sample_key = base_name.replace("_fastqc.zip", "")
        parts = sample_key.split("_")
        # Heuristic: if second token looks like read direction, keep it; otherwise treat whole as sample
        if len(parts) >= 2 and parts[1] in {"R1","R2","1","2"}:
            sample_name = parts[0]
            read_direction = parts[1] if parts[1] in {"R1","R2"} else ("R1" if parts[1]=="1" else "R2")
        else:
            sample_name = parts[0]
            read_direction = "R1"

        extracted_dir = extract_fastqc_zip(fastqc_zip)
        summary_file = find_summary_file(extracted_dir)
        
        if not summary_file:
            print(f"Warning: summary.txt not found in {fastqc_zip}. Skipping this sample.")
            continue
        
        status = parse_fastqc_summary(summary_file)
        sample_statuses.setdefault(sample_name, {})[read_direction] = status
        all_metrics.update(status.keys())
    
    # Convert set of all metrics to a sorted list for consistent ordering
    all_metrics = sorted(list(all_metrics))

    # Abbreviation table
    abbreviation_table = {
        'AC': 'Adapter Content',
        'BQ': 'Base Quality',
        'BGC': 'Per Base GC Content',
        'SGC': 'Per Sequence GC Content',
        'SN': 'Per Base N Content',
        'SL': 'Sequence Length Distribution',
        'SD': 'Sequence Duplication Levels',
        'Stat': 'Basic Statistics',
        'Kmer': 'Kmer Content',
        'BSC': 'Per Base Sequence Content',
        'UNK': 'Unknown Metric',
    }

    # Write the data to the main TSV file
    with open(output_tsv, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        # Write header row
        header = ["Sample"] + [f"{metric}_{read}" for metric in all_metrics for read in ["R1", "R2"]]
        writer.writerow(header)
        # Write data rows
        for sample, reads in sample_statuses.items():
            row = [sample]
            for metric in all_metrics:
                value_r1 = reads.get("R1", {}).get(metric, "")
                value_r2 = reads.get("R2", {}).get(metric, "")
                row.append(value_r1 if value_r1 != "N/A" else "")
                row.append(value_r2 if value_r2 != "N/A" else "")
            writer.writerow(row)
    
    # Write the abbreviation table to a separate legend file
    with open(legend_tsv, 'w', newline='') as legendfile:
        writer = csv.writer(legendfile, delimiter='\t')
        writer.writerow(["Abbreviation", "Description"])
        for code, description in abbreviation_table.items():
            writer.writerow([code, description])

def main():
    if len(sys.argv) != 4:
        print("Usage: python check_fastqc.py <output_tsv> <fastqc_directory> <legend_tsv>")
        sys.exit(1)
    
    output_tsv = sys.argv[1]
    fastqc_dir = sys.argv[2]
    legend_tsv = sys.argv[3]

    # Get a list of all FastQC .zip files in the specified directory
    fastqc_zips = [os.path.join(fastqc_dir, file) for file in os.listdir(fastqc_dir) if file.endswith('_fastqc.zip')]
    
    if not fastqc_zips:
        print(f"No FastQC .zip files found in {fastqc_dir}.")
        sys.exit(1)
    
    generate_report(fastqc_zips, output_tsv, legend_tsv)

if __name__ == "__main__":
    main()
