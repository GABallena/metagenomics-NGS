# Snakefile for Mobile Genetic Elements (MGE) analysis

# Configuration with defaults
# -------------------------
configfile: "configs/config.yaml"

# Default directory structure
config.setdefault("directories", {})
dirs = config["directories"]
dirs.setdefault("cleaned_reads", "data/cleaned_reads")
dirs.setdefault("merged_reads", "data/merged_reads") 
dirs.setdefault("output", "results")
dirs.setdefault("logs", "logs")
dirs.setdefault("plasmid_prediction", "results/plasmid_prediction")
dirs.setdefault("plasmid_assembly", "results/plasmid_assembly")
dirs.setdefault("mob_suite_output", "results/mob_suite")
dirs.setdefault("plasmidspades_output", "results/plasmidspades")
dirs.setdefault("gc_skew_analysis", "results/gc_skew")
dirs.setdefault("read_mapping", "results/read_mapping")

# Default conda environments
config.setdefault("envs", {})
envs = config["envs"]
envs.setdefault("plasmidfinder", "envs/plasmidfinder.yaml")
envs.setdefault("oritfinder", "envs/oritfinder.yaml")
envs.setdefault("tnppred", "envs/tnppred.yaml")
envs.setdefault("mob_suite", "envs/mob_suite.yaml")
envs.setdefault("plasmidspades", "envs/plasmidspades.yaml")
envs.setdefault("gapfiller", "envs/gapfiller.yaml")
envs.setdefault("recycler", "envs/recycler.yaml")

# Resource defaults
THREADS = config.get("threads", 4)
MEMORY = config.get("memory_mb", 8000)

# Create output directories
for dir in dirs.values():
    os.makedirs(dir, exist_ok=True)

# Define final output targets
rule all:
    input:
        mge_annotations = f"{dirs['output']}/mge_annotation/final_mges.gff",
        mge_types = f"{dirs['output']}/mge_classification/mge_types.tsv"


# Rule to identify MGEs using tools like RepeatMasker
rule identify_mges:
    input:
        genome = "data/genome.fasta"
    output:
        gff = "results/mge_annotation/raw_mges.gff"
    threads: 8
    conda:
        "envs/repeatmasker.yaml"
    shell:
        """
        RepeatMasker -pa {threads} -gff -species bacteria {input.genome} \
            > {output.gff}
        """

# Rule to classify MGEs
rule classify_mges:
    input:
        mges = "results/mge_annotation/raw_mges.gff"
    output:
        classified = "results/mge_annotation/final_mges.gff",
        types = "results/mge_classification/mge_types.tsv"
    conda:
        "envs/mge_tools.yaml"
    script:
        "scripts/classify_mges.py"

# Main workflow begins
rule all:
    input:
        f"{OUTPUT_DIR}/final_annotated_plasmids.fasta",
        f"{MOB_SUITE_OUTPUT_DIR}/mob_suite_summary.txt",
        f"{GC_SKEW_ANALYSIS_DIR}/gc_skew_plot.png"

# Top-down and bottom-up workflow (earlier steps) happen here, and this would link to cleaned reads from Kraken2 etc.

# Plasmid Prediction Step
rule predict_plasmids:
    input:
        cleaned_reads=f"{CLEANED_READS_DIR}/cleaned_reads.fastq.gz"
    output:
        predicted_plasmids=f"{PLASMID_PREDICTION_DIR}/predicted_plasmids.fasta"
    conda:
        RECYCLER_CONDA_ENV
    log:
        f"{LOG_DIR}/predict_plasmids.log"
    shell:
        """
        recycler -i {input.cleaned_reads} -o {output.predicted_plasmids} --threads {resources.threads} 2> {log}
        """

# Check with PlasmidFinder
rule check_plasmidfinder:
    input:
        predicted_plasmids=f"{PLASMID_PREDICTION_DIR}/predicted_plasmids.fasta"
    output:
        plasmidfinder_report=f"{PLASMID_PREDICTION_DIR}/plasmidfinder_report.txt"
    conda:
        PLASMIDFINDER_CONDA_ENV
    log:
        f"{LOG_DIR}/check_plasmidfinder.log"
    shell:
        """
        plasmidfinder.py -i {input.predicted_plasmids} -o {output.plasmidfinder_report} --threads {resources.threads} 2> {log}
        """

# OriTfinder Analysis
rule check_oritfinder:
    input:
        predicted_plasmids=f"{PLASMID_PREDICTION_DIR}/predicted_plasmids.fasta"
    output:
        oritfinder_report=f"{PLASMID_PREDICTION_DIR}/oritfinder_report.txt"
    conda:
        ORITFINDER_CONDA_ENV
    log:
        f"{LOG_DIR}/check_oritfinder.log"
    shell:
        """
        oritfinder.py -i {input.predicted_plasmids} -o {output.oritfinder_report} --threads {resources.threads} 2> {log}
        """

# MOB-Suite Analysis
rule run_mob_suite:
    input:
        predicted_plasmids=f"{PLASMID_PREDICTION_DIR}/predicted_plasmids.fasta"
    output:
        mob_suite_summary=f"{MOB_SUITE_OUTPUT_DIR}/mob_suite_summary.txt"
    conda:
        MOB_SUITE_CONDA_ENV
    log:
        f"{LOG_DIR}/run_mob_suite.log"
    shell:
        """
        mob_suite.py -i {input.predicted_plasmids} -o {output.mob_suite_summary} --threads {resources.threads} 2> {log}
        """

# PlasmidSPAdes Assembly
rule plasmid_assembly:
    input:
        predicted_plasmids=f"{PLASMID_PREDICTION_DIR}/predicted_plasmids.fasta"
    output:
        plasmid_assembly=f"{PLASMIDSPADES_OUTPUT_DIR}/plasmid_assembly.fasta"
    conda:
        PLASMIDSPADES_CONDA_ENV
    log:
        f"{LOG_DIR}/plasmid_assembly.log"
    shell:
        """
        plasmidspades.py -s {input.predicted_plasmids} -o {output.plasmid_assembly} --threads {resources.threads} 2> {log}
        """

# GapFiller to Close Gaps in Assembled Plasmids
rule close_gaps:
    input:
        plasmid_assembly=f"{PLASMIDSPADES_OUTPUT_DIR}/plasmid_assembly.fasta",
        reads=f"{MERGED_READS_DIR}/merged_reads.fastq.gz"
    output:
        final_plasmids=f"{OUTPUT_DIR}/final_annotated_plasmids.fasta"
    conda:
        GAPFILLER_CONDA_ENV
    log:
        f"{LOG_DIR}/close_gaps.log"
    shell:
        """
        GapFiller.pl -s {input.plasmid_assembly} -l libraries.txt -m 20 -o {output.final_plasmids} 2> {log}
        """

# GC Skew Analysis
rule calculate_gc_skew:
    input:
        final_plasmids=f"{OUTPUT_DIR}/final_annotated_plasmids.fasta"
    output:
        gc_skew_plot=f"{GC_SKEW_ANALYSIS_DIR}/gc_skew_plot.png"
    log:
        f"{LOG_DIR}/calculate_gc_skew.log"
    script:
        "scripts/gc_skew_analysis.py"


# Clustering of reads based on GC skew, tetranucleotide freq etc. 




# Normalization to 16S rRNA (FPKM)
rule normalize_to_16s:
    input:
        final_plasmids=f"{OUTPUT_DIR}/final_annotated_plasmids.fasta"
    output:
        fpkm_normalized=f"{OUTPUT_DIR}/fpkm_normalized_plasmids.txt"
    log:
        f"{LOG_DIR}/normalize_to_16s.log"
    script:
        "scripts/normalize_to_16s.py"




# Map Reads to Plasmids
rule map_reads_to_plasmids:
    input:
        reads=f"{MERGED_READS_DIR}/merged_reads.fastq.gz",
        final_plasmids=f"{OUTPUT_DIR}/final_annotated_plasmids.fasta"
    output:
        mapped_bam=f"{READ_MAPPING_DIR}/reads_mapped_to_plasmids.bam"
    conda:
        "path/to/bwa_env.yaml"
    log:
        f"{LOG_DIR}/map_reads_to_plasmids.log"
    shell:
        """
        bwa mem -t {resources.threads} {input.final_plasmids} {input.reads} | samtools view -bS - > {output.mapped_bam} 2> {log}
        samtools sort {output.mapped_bam} -o {output.mapped_bam}.sorted
        samtools index {output.mapped_bam}.sorted
        """

# Count Plasmid Reads
rule count_plasmid_reads:
    input:
        mapped_bam=f"{READ_MAPPING_DIR}/reads_mapped_to_plasmids.bam.sorted"
    output:
        plasmid_read_count=f"{READ_MAPPING_DIR}/plasmid_read_count.txt"
    log:
        f"{LOG_DIR}/count_plasmid_reads.log"
    shell:
        """
        samtools view -c {input.mapped_bam} > {output.plasmid_read_count}
        """

# Count Total Reads in Metagenomic Dataset
rule count_total_reads:
    input:
        reads=f"{MERGED_READS_DIR}/merged_reads.fastq.gz"
    output:
        total_read_count=f"{READ_MAPPING_DIR}/total_read_count.txt"
    log:
        f"{LOG_DIR}/count_total_reads.log"
    shell:
        """
        zcat {input.reads} | echo $((`wc -l`/4)) > {output.total_read_count}
        """

# Calculate Percentage of Plasmid Reads
rule calculate_plasmid_percentage:
    input:
        plasmid_read_count=f"{READ_MAPPING_DIR}/plasmid_read_count.txt",
        total_read_count=f"{READ_MAPPING_DIR}/total_read_count.txt"
    output:
        plasmid_percentage=f"{OUTPUT_DIR}/plasmid_percentage.txt"
    log:
        f"{LOG_DIR}/calculate_plasmid_percentage.log"
    script:
        "scripts/calculate_plasmid_percentage.py"
