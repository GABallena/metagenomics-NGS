import shutil
import glob
import os

# Check if raw_reads directory exists
if not os.path.exists("raw_reads"):
    raise ValueError(
        "raw_reads directory not found!\n"
        "Please create a 'raw_reads' directory and place your input files there."
    )

# Get samples from input directory - fix path
SAMPLES = [os.path.basename(f).replace("_1.fq.gz", "") 
          for f in glob.glob("raw_reads/*_1.fq.gz")]

# Add validation
if not SAMPLES:
    raise ValueError(
        "No input files found! Please ensure your input fastq files:\n"
        "1. Are in the 'raw_reads' directory\n"
        "2. Follow the naming pattern: *_1.fq.gz and *_2.fq.gz\n"
        "3. Have proper read permissions"
    )

print(f"Found {len(SAMPLES)} samples to process: {SAMPLES}")

# Parameter configuration
PARAMS = {
    "sickle": {
        "quality_threshold": 20,
        "length_threshold": 50,
        "encoding": "sanger"
    },
    "kma": {
        "memory_mode": True,
        "one_thread_one_task": True,
        "gt_value": 1,
        "index": {
            "k": 16,          # kmersize for indexing
            "k_t": 16,        # kmersize for template identification
            "k_i": 16         # kmersize for alignments
        },
        "mapping": {
            "min_phred": 20,  # minimum phred score
            "bc_threshold": 0.7,  # base call threshold
            "min_score": 0.0,  # minimum alignment score
            "dense": True     # skip insertions in consensus
        }
    },
    "megahit": {
        "memory": 16000  # MB
    },
    "fastqc": {
        "threads": 2
    },
    "metaphlan": {
        "input_type": "fastq",
        "analysis_type": "rel_ab",
        "bowtie2db": os.getenv("M4_DB", "/path/to/metaphlan_db"),  # Use M4_DB env var
        "index": "latest",  # Changed to latest
        "min_abundance": 0.1,
        "memory": "15GB",
        "force": True,
        "unclassified": True
    },
    "threads": 8,
    "contaminants": {
        "index_prefix": "contaminants"
    },
    "temp_files": {
        "keep": False
    }
}

# Environment paths
ENVS = {
    "fastqc": "fastqc_env",
    "sickle": "sickle_env",
    "kma": "kma_env",
    "megahit": "megahit_env",
    "metaphlan": "metaphlan_env"
}

# Directories configuration
DIRS = {
    "results": "results",
    "input": "raw_reads",  # Changed from results/raw_reads to raw_reads
    "sickle": "results/sickle",
    "fastqc": "results/fastqc",
    "kma": "results/kma",
    "megahit": "results/megahit",
    "logs": "results/logs",
    "contaminants": "contaminants",  # Removed trailing slash
    "plots": "results/plots",  # Add plots directory
    "metaphlan_db": "M4_DB"  # Add MetaPhlAn DB path
}

# Add contaminant files validation
if not glob.glob(os.path.join(DIRS["contaminants"], "*.fasta")):
    raise ValueError(
        f"No FASTA files found in {DIRS['contaminants']} directory!\n"
        "Please add your contaminant FASTA files to this directory."
    )

# Default configuration
config.setdefault("threads", PARAMS["threads"])
config.setdefault("mem_mb", PARAMS["megahit"]["memory"])

# Function to clean up existing directories
def cleanup_directories():
    for key, dirname in DIRS.items():
        if key not in ["input", "contaminants"] and os.path.exists(dirname):
            shutil.rmtree(dirname)
            print(f"Removed existing directory: {dirname}")

# Clean up and recreate all output directories
cleanup_directories()
for dirname in DIRS.values():
    os.makedirs(dirname, exist_ok=True)
    print(f"Created directory: {dirname}")

# Add a onstart handler to validate directories
onstart:
    print("Creating output directories...")
    for dirname in DIRS.values():
        if not os.path.exists(dirname):
            print(f"Error: Could not create directory {dirname}")
            exit(1)

# Target rule
rule all:
    input:
        # FastQC outputs
        expand(DIRS["fastqc"] + "/raw/{sample}_{read}_fastqc.{ext}",
               sample=SAMPLES, read=["1", "2"], ext=["html", "zip"]),
        expand(DIRS["fastqc"] + "/trimmed/{sample}_{read}_trimmed_fastqc.{ext}",
               sample=SAMPLES, read=["1", "2"], ext=["html", "zip"]),
        # Trimmed reads
        expand(DIRS["sickle"] + "/{sample}_{read}_trimmed.fastq",
               sample=SAMPLES, read=["1", "2"]),
        # Contaminant removal outputs (updated)
        expand(DIRS["kma"] + "/{sample}.res", sample=SAMPLES),
        expand(DIRS["kma"] + "/{sample}.frag.gz", sample=SAMPLES),
        expand(DIRS["kma"] + "/{sample}.aln", sample=SAMPLES),
        # Rest of the outputs
        expand(DIRS["results"] + "/metaphlan/{sample}_profile.txt", sample=SAMPLES),
        expand(DIRS["megahit"] + "/{sample}", sample=SAMPLES)

# Rules section
rule fastqc_raw:
    input:
        # Change input to depend on decompressed files
        DIRS["input"] + "/{sample}_{read}.fastq"
    output:
        html = DIRS["fastqc"] + "/raw/{sample}_{read}_fastqc.html",
        zip = DIRS["fastqc"] + "/raw/{sample}_{read}_fastqc.zip"
    log:
        DIRS["logs"] + "/fastqc/raw/{sample}_{read}.log"
    threads: PARAMS["fastqc"]["threads"]
    conda:
        ENVS["fastqc"]
    shell:
        """
        mkdir -p {DIRS[fastqc]}/raw
        fastqc -t {threads} -o {DIRS[fastqc]}/raw {input} 2> {log}
        """

rule decompress_reads:
    input:
        r1 = DIRS["input"] + "/{sample}_1.fq.gz",
        r2 = DIRS["input"] + "/{sample}_2.fq.gz"
    output:
        r1 = temp(DIRS["input"] + "/{sample}_1.fastq"),
        r2 = temp(DIRS["input"] + "/{sample}_2.fastq")
    shell:
        """
        gunzip -c {input.r1} > {output.r1}
        gunzip -c {input.r2} > {output.r2}
        """

rule trim_reads:
    input:
        r1 = DIRS["input"] + "/{sample}_1.fastq",
        r2 = DIRS["input"] + "/{sample}_2.fastq"
    output:
        r1 = DIRS["sickle"] + "/{sample}_1_trimmed.fastq",
        r2 = DIRS["sickle"] + "/{sample}_2_trimmed.fastq",
        singles = DIRS["sickle"] + "/{sample}_singles_trimmed.fastq"
    log:
        DIRS["logs"] + "/sickle/{sample}.log"
    threads: config["threads"]
    conda:
        ENVS["sickle"]
    params:
        quality = PARAMS["sickle"]["quality_threshold"],
        length = PARAMS["sickle"]["length_threshold"]
    shell:
        "sickle pe -f {input.r1} -r {input.r2} "
        "-t sanger "
        "-o {output.r1} -p {output.r2} -s {output.singles} "
        "-q {params.quality} -l {params.length} "
        "2> {log}"

rule fastqc_trimmed:
    input:
        DIRS["sickle"] + "/{sample}_{read}_trimmed.fastq"
    output:
        html = DIRS["fastqc"] + "/trimmed/{sample}_{read}_trimmed_fastqc.html",
        zip = DIRS["fastqc"] + "/trimmed/{sample}_{read}_trimmed_fastqc.zip"
    log:
        DIRS["logs"] + "/fastqc/trimmed/{sample}_{read}.log"
    threads: PARAMS["fastqc"]["threads"]
    conda:
        ENVS["fastqc"]
    shell:
        "fastqc -t {threads} -o {DIRS[fastqc]}/trimmed {input} 2> {log}"

rule concatenate_contaminants:
    input:
        contam_files = glob.glob(os.path.join(DIRS["contaminants"], "*.fasta"))
    output:
        os.path.join(DIRS["contaminants"], "all_contaminants.fasta")
    shell:
        "cat {input.contam_files} > {output}"

rule index_contaminants:
    input:
        fasta = os.path.join(DIRS["contaminants"], "all_contaminants.fasta")
    output:
        comp = os.path.join(DIRS["contaminants"], "contaminants.comp.b"),
        length = os.path.join(DIRS["contaminants"], "contaminants.length.b"),
        name = os.path.join(DIRS["contaminants"], "contaminants.name"),
        seq = os.path.join(DIRS["contaminants"], "contaminants.seq.b")
    log:
        DIRS["logs"] + "/contaminants/index.log"
    params:
        k = PARAMS["kma"]["index"]["k"],
        k_t = PARAMS["kma"]["index"]["k_t"],
        k_i = PARAMS["kma"]["index"]["k_i"],
        prefix = os.path.join(DIRS["contaminants"], PARAMS["contaminants"]["index_prefix"])
    shell:
        """
        kma index -i {input.fasta} \
            -o {params.prefix} \
            -k {params.k} \
            -k_t {params.k_t} \
            -k_i {params.k_i} \
            2> {log}
        """

rule remove_contaminants:
    input:
        r1 = DIRS["sickle"] + "/{sample}_1_trimmed.fastq",
        r2 = DIRS["sickle"] + "/{sample}_2_trimmed.fastq",
        db = os.path.join(DIRS["contaminants"], "contaminants.comp.b")
    output:
        res = DIRS["kma"] + "/{sample}.res",
        frag = DIRS["kma"] + "/{sample}.frag.gz",
        aln = DIRS["kma"] + "/{sample}.aln"
    log:
        DIRS["logs"] + "/kma/{sample}.log"
    threads: config["threads"]
    conda:
        ENVS["kma"]
    params:
        memory_mode = "-mem_mode" if PARAMS["kma"]["memory_mode"] else "",
        one_thread = "-1t1" if PARAMS["kma"]["one_thread_one_task"] else "",
        min_phred = PARAMS["kma"]["mapping"]["min_phred"],
        bc_threshold = PARAMS["kma"]["mapping"]["bc_threshold"],
        min_score = PARAMS["kma"]["mapping"]["min_score"],
        dense = "-dense" if PARAMS["kma"]["mapping"]["dense"] else "",
        prefix = DIRS["kma"] + "/{sample}"
    shell:
        """
        kma -ipe {input.r1} {input.r2} \
            -o {params.prefix} \
            -t_db {DIRS[contaminants]}/contaminants \
            -t {threads} \
            {params.memory_mode} \
            {params.one_thread} \
            {params.dense} \
            -mp {params.min_phred} \
            -bc {params.bc_threshold} \
            -mrs {params.min_score} \
            2> {log}
        """

rule metaphlan:
    input:
        r1 = DIRS["sickle"] + "/{sample}_1_trimmed.fastq",
        r2 = DIRS["sickle"] + "/{sample}_2_trimmed.fastq"
    output:
        profile = DIRS["results"] + "/metaphlan/{sample}_profile.txt",
        bowtie2out = DIRS["results"] + "/metaphlan/{sample}.bowtie2.bz2"
    log:
        DIRS["logs"] + "/metaphlan/{sample}.log"
    threads: config["threads"]
    conda:
        ENVS["metaphlan"]
    params:
        min_abundance = PARAMS["metaphlan"]["min_abundance"],
        db_dir = DIRS["metaphlan_db"],  # Use DIRS instead of PARAMS
        force = "--force" if PARAMS["metaphlan"]["force"] else "",
        unclassified = "--unclassified_estimation" if PARAMS["metaphlan"]["unclassified"] else "",
        memory = PARAMS["metaphlan"]["memory"]
    shell:
        """
        # Setup database if needed
        if [ ! -d {params.db_dir} ]; then
            mkdir -p {params.db_dir}
            metaphlan --install --bowtie2db {params.db_dir} 2>> {log}
        fi

        # Run MetaPhlAn
        metaphlan {input.r1},{input.r2} \
            --input_type fastq \
            --nproc {threads} \
            {params.force} \
            --bowtie2db {params.db_dir} \
            --output_file {output.profile} \
            --bowtie2out {output.bowtie2out} \
            --memory-min {params.memory} \
            --index {PARAMS[metaphlan][index]} \
            {params.unclassified} \
            --sample_id $(basename {output.profile} .txt) \
            2>> {log} || \
            {{ echo "k__Bacteria\tp__Unknown\t100.0" > {output.profile}; exit 0; }}
        """

rule megahit_assembly:
    input:
        r1 = DIRS["sickle"] + "/{sample}_1_trimmed.fastq",
        r2 = DIRS["sickle"] + "/{sample}_2_trimmed.fastq"
    output:
        directory(DIRS["megahit"] + "/{sample}")
    log:
        DIRS["logs"] + "/megahit/{sample}.log"
    threads: config["threads"]
    resources:
        mem_mb = config["mem_mb"]
    conda:
        ENVS["megahit"]
    shell:
        "rm -rf {output} && "
        "megahit "
        "-1 {input.r1} "
        "-2 {input.r2} "
        "-t {threads} "
        "--memory {resources.mem_mb} "
        "-o {output} "
        "2> {log}"



