import os
import glob
from snakemake.utils import min_version

# Check Snakemake version
min_version("6.0")

# Constants
READS_DIR = "trimmed_reads"
FASTUNIQ_DIR = "processed_reads"
ASSEMBLY_DIR = "spades_output"
LONG_READS_DIR = "long_reads"
LOG_DIR = "logs"
BENCHMARK_DIR = "benchmarks"
FASTUNIQ_ENV = "fastuniq_env"
SPADES_ENV = "spades_env"
PANDASEQ_ENV = "pandaseq_env"

# Resource configurations
FASTUNIQ_THREADS = 4
SPADES_THREADS = 8
MAX_MEMORY_GB = 16
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
for dir in [FASTUNIQ_DIR, ASSEMBLY_DIR, LONG_READS_DIR, LOG_DIR, BENCHMARK_DIR]:
    os.makedirs(dir, exist_ok=True)

SAMPLES = get_samples()

rule all:
    input:
        expand(f"{FASTUNIQ_DIR}/{{sample}}_{{read}}.fastq", sample=SAMPLES, read=[1, 2]),
        expand(f"{ASSEMBLY_DIR}/{{sample}}_assembly/{{sample}}_final_filtered_contigs.fa", sample=SAMPLES),
        expand(f"{LONG_READS_DIR}/{{sample}}_longreads.fastq", sample=SAMPLES),
        expand(f"{LONG_READS_DIR}/{{sample}}_quality_metrics.txt", sample=SAMPLES)

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

# Rule to deduplicate reads with FastUniq
rule deduplicate_reads:
    input:
        read1_unzipped=f"{READS_DIR}/{{sample}}_R1_paired.fastq",
        read2_unzipped=f"{READS_DIR}/{{sample}}_R2_paired.fastq"
    output:
        read1_clean=f"{FASTUNIQ_DIR}/{{sample}}_1.fastq",
        read2_clean=f"{FASTUNIQ_DIR}/{{sample}}_2.fastq"
    threads: FASTUNIQ_THREADS
    conda: FASTUNIQ_ENV
    shell:
        """
        echo -e "{input.read1_unzipped}\\n{input.read2_unzipped}" > {FASTUNIQ_DIR}/file_list.txt
        fastuniq -i {FASTUNIQ_DIR}/file_list.txt -t q -o {output.read1_clean} -p {output.read2_clean}
        """

rule spades:
    input:
        read1=f"{FASTUNIQ_DIR}/{{sample}}_1.fastq",
        read2=f"{FASTUNIQ_DIR}/{{sample}}_2.fastq",
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
    conda: SPADES_ENV
    benchmark:
        f"{BENCHMARK_DIR}/spades_{{sample}}.txt"
    log:
        main=f"{LOG_DIR}/spades_{{sample}}.log",
        stats=f"{LOG_DIR}/spades_stats_{{sample}}.log"
    shell:
        """
        chmod +x {input.script}
        {input.script} \
            {input.read1} \
            {input.read2} \
            {params.outdir} \
            {output.filtered_contigs} \
            {threads} \
            {log.main} \
            {params.min_len} || {{ echo "SPAdes assembly failed" >&2; exit 1; }}
            
        # Generate stats if assembly succeeded
        seqkit stats {output.filtered_contigs} > {output.stats} 2> {log.stats}
        """

# Rule: Merge Reads with PANDAseq
rule run_pandaseq:
    input:
        read1=f"{FASTUNIQ_DIR}/{{sample}}_1.fastq",
        read2=f"{FASTUNIQ_DIR}/{{sample}}_2.fastq",
        script="scripts/pandaseq.sh"
    output:
        merged_reads=f"{LONG_READS_DIR}/{{sample}}_longreads.fastq",
        quality_metrics=f"{LONG_READS_DIR}/{{sample}}_quality_metrics.txt"
    params:
        min_overlap=10,
        quality_threshold=0.9,
        outdir=lambda wildcards, output: os.path.dirname(output.merged_reads)
    conda: PANDASEQ_ENV
    threads: 4
    resources:
        mem_mb=4000,
        runtime=60
    benchmark:
        f"{BENCHMARK_DIR}/pandaseq_{{sample}}.txt"
    log:
        main=f"{LOG_DIR}/pandaseq_{{sample}}.log",
        stdout=f"{LOG_DIR}/pandaseq_{{sample}}_stdout.log",
        stderr=f"{LOG_DIR}/pandaseq_{{sample}}_stderr.log"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log.main})
        chmod +x {input.script}
        
        SAMPLE={wildcards.sample} \\
        READ1={input.read1} \\
        READ2={input.read2} \\
        MERGED_READS={output.merged_reads} \\
        QUALITY_METRICS={output.quality_metrics} \\
        LOG={log.main} \\
        MIN_OVERLAP={params.min_overlap} \\
        QUALITY_THRESHOLD={params.quality_threshold} \\
        bash {input.script} > {log.stdout} 2> {log.stderr}
        """



