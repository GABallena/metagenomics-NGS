import os
import glob
from snakemake.utils import min_version

# Check Snakemake version
min_version("6.0")

# Constants
READS_DIR = "trimmed_reads"
ASSEMBLY_DIR = "spades_output"
LOG_DIR = "logs"
BENCHMARK_DIR = "benchmarks"
SPADES_ENV = "spades_env"

# Resource configurations
SPADES_THREADS = 8
MAX_MEMORY_GB = 24
MIN_CONTIG_LEN = 200

# Input validation function
def validate_input_files():
    files = glob.glob(f"{READS_DIR}/*_R1_paired.fastq.gz")
    if not files:
        raise ValueError(f"No input files found in {READS_DIR}")
    return sorted(files)

# Get samples with validation
def get_samples():
    files = validate_input_files()
    return [os.path.basename(f).replace('_R1_paired.fastq.gz', '') for f in files]

# Create directories
for dir in [ASSEMBLY_DIR, LOG_DIR, BENCHMARK_DIR]:
    os.makedirs(dir, exist_ok=True)

SAMPLES = get_samples()

rule all:
    input:
        expand(f"{ASSEMBLY_DIR}/{{sample}}_assembly/{{sample}}_final_filtered_contigs.fa", sample=SAMPLES)

rule unzip_fastq:
    input:
        read1=f"{READS_DIR}/{{sample}}_R1_paired.fastq.gz",
        read2=f"{READS_DIR}/{{sample}}_R2_paired.fastq.gz"
    output:
        read1_unzipped=temp(f"{READS_DIR}/{{sample}}_R1_paired.fastq"),
        read2_unzipped=temp(f"{READS_DIR}/{{sample}}_R2_paired.fastq")
    benchmark:
        f"{BENCHMARK_DIR}/unzip_{{sample}}.txt"
    log:
        f"{LOG_DIR}/unzip_{{sample}}.log"
    shell:
        """
        gunzip -c {input.read1} > {output.read1_unzipped} 2> {log}
        gunzip -c {input.read2} > {output.read2_unzipped} 2>> {log}
        """

rule spades:
    input:
        read1=f"{READS_DIR}/{{sample}}_R1_paired.fastq.gz",
        read2=f"{READS_DIR}/{{sample}}_R2_paired.fastq.gz",
        script="scripts/spades_binning.sh"
    output:
        filtered_contigs=f"{ASSEMBLY_DIR}/{{sample}}_assembly/{{sample}}_final_filtered_contigs.fa",
        assembly_graph=f"{ASSEMBLY_DIR}/{{sample}}_assembly/assembly_graph.fastg",
        stats=f"{ASSEMBLY_DIR}/{{sample}}_assembly/assembly_stats.txt"
    params:
        outdir=f"{ASSEMBLY_DIR}/{{sample}}_assembly",
        min_len=MIN_CONTIG_LEN
    threads: SPADES_THREADS
    resources:
        mem_mb=lambda wildcards, threads: threads * 4000
    conda: "spades_env" 
    benchmark:
        f"{BENCHMARK_DIR}/spades_{{sample}}.txt"
    log:
        main=f"{LOG_DIR}/spades_{{sample}}.log",
        stats=f"{LOG_DIR}/spades_stats_{{sample}}.log"
    shell:
        """
        # Ensure SPAdes is available
        command -v spades.py >/dev/null 2>&1 || {{ echo "SPAdes is not installed"; exit 1; }}
        
        chmod +x {input.script}
        {input.script} \\
            {input.read1} \\
            {input.read2} \\
            {params.outdir} \\
            {output.filtered_contigs} \\
            {threads} \\
            {log.main} \\
            {params.min_len} || {{ echo "SPAdes assembly failed" >&2; exit 1; }}
            
        # Generate stats if assembly succeeded
        seqkit stats {output.filtered_contigs} > {output.stats} 2> {log.stats}
        """



