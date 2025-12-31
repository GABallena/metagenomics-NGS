configfile: "config.yaml"
import os
from pathlib import Path

# Configuration dictionary
config = {
    "dirs": {
        "raw_reads": os.path.dirname(config.get("r1", "")),
        "trimmed": "trimmed_reads",
        "fastqc": "fastqc_output",
        "kraken": "Kraken/kraken_output",
        "bracken": "Kraken/bracken_output",
        "results": "aggregated_results",
        "logs": "logs"
    },
    "params": {
        "trim": {
            "adapters": "adapters.fasta",
            "slidingwindow": "4:15",
            "leading": 3,
            "trailing": 3,
            "minlen": 36
        }
    },
    "threads": {
        "trim": 16,
        "fastqc": 8,
        "kraken": 16,
        "diversity": 16
    }
}

# Now safe to define sample
raw_reads = config.get("r1", "").replace("_1.fq.gz", "").replace("_2.fq.gz", "")
sample_name = os.path.basename(raw_reads)

# Validate input files
def validate_inputs():
    required_files = [config.get("r1"), config.get("r2"), "adapters.fasta"]
    for f in required_files:
        if not f or not Path(f).is_file():
            raise ValueError(f"Required file {f} not found")

# Create necessary directories
for dir in config["dirs"].values():
    Path(dir).mkdir(parents=True, exist_ok=True)

# Validate inputs
validate_inputs()

# Single sample processing
SAMPLES = [sample_name]

def get_mem_mb(wildcards, attempt):
    return 4000 * attempt

rule all:
    input:
        expand(f"{config['dirs']['fastqc']}/{{sample}}_R1_paired_fastqc.zip", sample=SAMPLES),
        expand(f"{config['dirs']['fastqc']}/{{sample}}_R2_paired_fastqc.zip", sample=SAMPLES),
        expand(f"{config['dirs']['fastqc']}/{{sample}}_R1_paired_fastqc.html", sample=SAMPLES),
        expand(f"{config['dirs']['fastqc']}/{{sample}}_R2_paired_fastqc.html", sample=SAMPLES),
        expand(f"{config['dirs']['trimmed']}/{{sample}}_R1_paired.fq.gz", sample=SAMPLES),
        expand(f"{config['dirs']['trimmed']}/{{sample}}_R2_paired.fq.gz", sample=SAMPLES),
        expand(f"{config['dirs']['kraken']}/{{sample}}.k2report", sample=SAMPLES),
        expand(f"{config['dirs']['kraken']}/{{sample}}.kraken2", sample=SAMPLES),
        expand(f"{config['dirs']['bracken']}/{{sample}}.bracken", sample=SAMPLES),
        expand(f"{config['dirs']['bracken']}/{{sample}}.breport", sample=SAMPLES)

rule trimming:
    input:
        R1=lambda wildcards: f"{config['dirs']['raw_reads']}/{wildcards.sample}_R1_001.fastq.gz",
        R2=lambda wildcards: f"{config['dirs']['raw_reads']}/{wildcards.sample}_R2_001.fastq.gz"
    output:
        R1_trimmed=f"{config['dirs']['trimmed']}/{{sample}}_R1_paired.fq.gz",
        R2_trimmed=f"{config['dirs']['trimmed']}/{{sample}}_R2_paired.fq.gz",
        R1_unpaired=temp(f"{config['dirs']['trimmed']}/{{sample}}_R1_unpaired.fq.gz"),
        R2_unpaired=temp(f"{config['dirs']['trimmed']}/{{sample}}_R2_unpaired.fq.gz")
    threads: config["threads"]["trim"]
    resources:
        mem_mb=get_mem_mb
    conda: "trimmomatic_env"
    log:
        f"{config['dirs']['logs']}/trimming/{{sample}}.log"
    params:
        adapters=config["params"]["trim"]["adapters"],
        leading=config["params"]["trim"]["leading"],
        trailing=config["params"]["trim"]["trailing"],
        slidingwindow=config["params"]["trim"]["slidingwindow"],
        minlen=config["params"]["trim"]["minlen"]
    shell:
        """
        mkdir -p {config['dirs']['trimmed']} {config['dirs']['logs']}/trimming
        trimmomatic PE -threads {threads} -phred33 {input.R1} {input.R2} \
        {output.R1_trimmed} {output.R1_unpaired} \
        {output.R2_trimmed} {output.R2_unpaired} \
        ILLUMINACLIP:{params.adapters}:2:30:10 \
        LEADING:{params.leading} TRAILING:{params.trailing} \
        SLIDINGWINDOW:{params.slidingwindow} MINLEN:{params.minlen} \
        2> {log}
        """

rule fastqc:
    input:
        R1=f"{config['dirs']['trimmed']}/{{sample}}_R1_paired.fq.gz",
        R2=f"{config['dirs']['trimmed']}/{{sample}}_R2_paired.fq.gz"
    output:
        fastqc_zip1=f"{config['dirs']['fastqc']}/{{sample}}_R1_paired_fastqc.zip",
        fastqc_zip2=f"{config['dirs']['fastqc']}/{{sample}}_R2_paired_fastqc.zip",
        fastqc_html1=f"{config['dirs']['fastqc']}/{{sample}}_R1_paired_fastqc.html",
        fastqc_html2=f"{config['dirs']['fastqc']}/{{sample}}_R2_paired_fastqc.html"
    threads: config["threads"]["fastqc"]
    resources:
        mem_mb=4000
    conda: "fastqc_env"
    log:
        f"{config['dirs']['logs']}/fastqc/{{sample}}.log"
    shell:
        """
        mkdir -p {config['dirs']['fastqc']} {config['dirs']['logs']}/fastqc
        fastqc {input.R1} {input.R2} -o {config['dirs']['fastqc']} 2> {log}
        """

rule kraken_bracken:
    input:
        R1=f"{config['dirs']['trimmed']}/{{sample}}_R1_paired.fq.gz",
        R2=f"{config['dirs']['trimmed']}/{{sample}}_R2_paired.fq.gz"
    output:
        kraken_report=f"{config['dirs']['kraken']}/{{sample}}.k2report",
        kraken_output=f"{config['dirs']['kraken']}/{{sample}}.kraken2",
        bracken_output=f"{config['dirs']['bracken']}/{{sample}}.bracken",
        bracken_report=f"{config['dirs']['bracken']}/{{sample}}.breport"
    threads: config["threads"]["kraken"]
    resources:
        mem_mb=8000
    conda: "taxonomy_env"
    log:
        f"{config['dirs']['logs']}/kraken_bracken/{{sample}}.log"
    shell:
        """
        mkdir -p {config['dirs']['kraken']} {config['dirs']['bracken']} {config['dirs']['logs']}/kraken_bracken
        kraken2 --db Kraken/k2_db --threads {threads} \
        --paired {input.R1} {input.R2} --report {output.kraken_report} > {output.kraken_output}
        
        bracken -d Kraken/k2_db -i {output.kraken_report} -o {output.bracken_output}
        
        cp {output.bracken_output} {output.bracken_report} 2> {log}
        """
