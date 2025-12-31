"""
Experimental script used to validate trimming and QC parameters
during early pipeline development.

Not intended for HPC-scale execution.
Superseded by Snakemake-based implementation.
"""
import os
import subprocess

# Directories and file paths
RAW_READS_DIR = "raw_reads"
TRIMMED_READS_DIR = "trimmed_reads"
FASTQC_REPORTS_DIR = "fastqc_reports"
TRIMMOMATIC_ADAPTERS = "adapters.fa"  # Ensure this file is in the same directory as the script

# Conda environments
CONDA_ENV_TRIMMOMATIC = "trimmomatic_env"
CONDA_ENV_FASTQC = "fastqc_env"

def run_command(cmd, conda_env=None):
    """
    Runs a shell command with optional Conda environment activation.
    """
    if conda_env:
        cmd = f"source ~/miniconda3/etc/profile.d/conda.sh && conda activate {conda_env} && {cmd}"
    subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')

def compress_with_seqtk(input_file, output_file):
    """
    Compresses a .fastq file to .fastq.gz using seqtk.
    """
    seqtk_cmd = f"seqtk seq {input_file} | gzip > {output_file}"
    run_command(seqtk_cmd)

def trim_reads(sample):
    """
    Trims paired-end reads using Trimmomatic and saves the results in the trimmed_reads directory.
    """
    read_1_path = f"{RAW_READS_DIR}/{sample}_R1_001.fastq.gz"
    read_2_path = f"{RAW_READS_DIR}/{sample}_R2_001.fastq.gz"
    paired_1 = f"{TRIMMED_READS_DIR}/{sample}_1_paired.fastq"
    unpaired_1 = f"{TRIMMED_READS_DIR}/{sample}_1_unpaired.fastq"
    paired_2 = f"{TRIMMED_READS_DIR}/{sample}_2_paired.fastq"
    unpaired_2 = f"{TRIMMED_READS_DIR}/{sample}_2_unpaired.fastq"

    if not os.path.isfile(read_1_path) or not os.path.isfile(read_2_path):
        print(f"Error: One or both files are missing for sample {sample}")
        return False

    os.makedirs(TRIMMED_READS_DIR, exist_ok=True)

    trimmomatic_cmd = (
        f"trimmomatic PE -threads 4 {read_1_path} {read_2_path} "
        f"{paired_1} {unpaired_1} {paired_2} {unpaired_2} "
        f"ILLUMINACLIP:{TRIMMOMATIC_ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 "
        "SLIDINGWINDOW:4:20 MINLEN:36"
    )
    run_command(trimmomatic_cmd, conda_env=CONDA_ENV_TRIMMOMATIC)
    return True

def run_fastqc(sample):
    """
    Runs FastQC on trimmed reads and stores the reports in the fastqc_reports directory.
    """
    paired_1 = f"{TRIMMED_READS_DIR}/{sample}_1_paired.fastq"
    paired_2 = f"{TRIMMED_READS_DIR}/{sample}_2_paired.fastq"

    if not os.path.isfile(paired_1) or not os.path.isfile(paired_2):
        print(f"Skipping FastQC for {sample}: Trimmed files not found.")
        return

    os.makedirs(FASTQC_REPORTS_DIR, exist_ok=True)

    fastqc_cmd = f"fastqc {paired_1} {paired_2} --outdir {FASTQC_REPORTS_DIR}"
    run_command(fastqc_cmd, conda_env=CONDA_ENV_FASTQC)

def compress_trimmed_reads(samples):
    """
    Compresses all paired and unpaired .fastq files in the trimmed_reads directory.
    """
    for sample in samples:
        paired_1 = f"{TRIMMED_READS_DIR}/{sample}_1_paired.fastq"
        paired_2 = f"{TRIMMED_READS_DIR}/{sample}_2_paired.fastq"
        unpaired_1 = f"{TRIMMED_READS_DIR}/{sample}_1_unpaired.fastq"
        unpaired_2 = f"{TRIMMED_READS_DIR}/{sample}_2_unpaired.fastq"

        for f in [paired_1, paired_2, unpaired_1, unpaired_2]:
            if os.path.isfile(f):
                compress_with_seqtk(f, f"{f}.gz")
                os.remove(f)

def generate_summary_report(samples):
    """
    Generates a summary report for processed samples.
    """
    summary_file = "cleaning_results/summary_report.txt"
    os.makedirs("cleaning_results", exist_ok=True)

    with open(summary_file, "w") as f:
        f.write("Sample\tStatus\n")
        for sample in samples:
            paired_1 = f"{TRIMMED_READS_DIR}/{sample}_1_paired.fastq.gz"
            paired_2 = f"{TRIMMED_READS_DIR}/{sample}_2_paired.fastq.gz"

            if os.path.isfile(paired_1) and os.path.isfile(paired_2):
                f.write(f"{sample}\tProcessed\n")
            else:
                f.write(f"{sample}\tFailed\n")

    print(f"Summary report generated: {summary_file}")

def process_samples(samples):
    """
    Processes each sample: trims reads, runs FastQC, and compresses the trimmed reads.
    """
    for sample in samples:
        print(f"Processing {sample}...")

        if not trim_reads(sample):
            print(f"Error trimming reads for {sample}. Skipping FastQC.")
            continue

        print(f"Trimming completed for {sample}.")
        run_fastqc(sample)
        print(f"FastQC completed for {sample}.")

    print("Compressing trimmed reads...")
    compress_trimmed_reads(samples)
    print("Compression completed.")
    generate_summary_report(samples)

if __name__ == "__main__":
    # Automatically derive samples from filenames
    samples = sorted(set(
        f.rsplit("_R", 1)[0] for f in os.listdir(RAW_READS_DIR) if f.endswith(".fastq.gz")
    ))
    process_samples(samples)
