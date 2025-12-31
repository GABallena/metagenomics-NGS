# Define global parameters
# ------------------------

# Directory paths
CLEANED_READS_DIR = "cleaned_reads"
MERGED_READS_DIR = "merged_reads"
CARD_DB_DIR = "databases/CARD"
OUTPUT_DIR = "output"
LOG_DIR = "logs"
PLASPREDICT_OUTPUT_DIR = "output/plaspredict"
PFAM_DB = "databases/Pfam/Pfam-A.hmm"

# Tool-specific resources
# PEAR
PEAR_THREADS = 8
PEAR_MEMORY = 16000
PEAR_CONDA_ENV = "envs/pear_env.yaml"

# KMA 
KMA_THREADS = 8
KMA_MEMORY = 32000
KMA_CONDA_ENV = "envs/kma_env.yaml"
KMA_ITERATIONS = 5

# PANDAseq
PANDASEQ_THREADS = 8
PANDASEQ_MEMORY = 16000
PANDASEQ_CONDA_ENV = "envs/pandaseq_env.yaml"

# metaSPAdes
SPADES_THREADS = 16
SPADES_MEMORY = 64000
SPADES_CONDA_ENV = "envs/spades_env.yaml"

# RGI
RGI_THREADS = 8
RGI_MEMORY = 32000
RGI_CONDA_ENV = "envs/rgi_env.yaml"

# Kraken2
KRAKEN2_THREADS = 16
KRAKEN2_MEMORY = 64000
KRAKEN2_CONDA_ENV = "envs/kraken2_env.yaml"

# PlasPredict
PLASPREDICT_THREADS = 8
PLASPREDICT_MEMORY = 32000
PLASPREDICT_CONDA_ENV = "envs/plaspredict_env.yaml"

# HMMER (for PFAM and TnpPred)
HMMER_THREADS = 8
HMMER_MEMORY = 16000
HMMER_CONDA_ENV = "envs/hmmer_env.yaml"

# Snakemake rules
# ---------------

rule all:
    input:
        f"{LOG_DIR}/final_kraken2_report.txt",
        f"{OUTPUT_DIR}/final_filtered_contigs.fasta",
        f"{PLASPREDICT_OUTPUT_DIR}/",
        f"{OUTPUT_DIR}/categorized_MGEs.txt"

# Forward Pipeline Rules

# Step 1: PEAR for merging paired-end reads
rule run_pear:
    input:
        r1="trimmed_reads/{sample}_R1_paired.fastq.gz",
        r2="trimmed_reads/{sample}_R2_paired.fastq.gz"
    output:
        merged="pair_merged/{sample}_merged.fastq.gz",
        discarded="pair_merged/{sample}_discarded.fastq.gz"
    log:
        "logs/{sample}_pear.log"
    conda:
        "envs/pear_env.yaml"
    shell:
        """
        pear -f {input.r1} -r {input.r2} -o pair_merged/{wildcards.sample}_ \
            2> {log}
        """

# Step 2: Convert FASTQ to FASTA using seqtk (for mate-pairs)
rule convert_fastq_to_fasta:
    input:
        fastq_gz="pair_merged/{sample}_merged.fastq.gz"
    output:
        fasta_gz="pair_merged/{sample}_merged.fasta.gz"
    log:
        "logs/{sample}_convert_fastq_to_fasta.log"
    conda:
        "envs/seqtk_env.yaml"
    shell:
        """
        seqtk seq -A {input.fastq_gz} > {output.fasta_gz} 2> {log}
        """

# Step 2.1: Translate CARD sequences into protein sequences
rule translate_card_sequences:
    input:
        card_fasta="path/to/card_sequences.fasta"
    output:
        translated_proteins="card_translated/{sample}_translated_proteins.fasta"
    log:
        "logs/{sample}_translate_card_sequences.log"
    conda:
        "envs/emboss.yaml"
    shell:
        """
        transeq -sequence {input.card_fasta} -outseq {output.translated_proteins} 2> {log}
        """

# Step 2.2: Reverse translate protein sequences back to nucleotide sequences
rule reverse_translate_proteins:
    input:
        translated_proteins="card_translated/{sample}_translated_proteins.fasta"
    output:
        reverse_translated_fasta="card_translated/{sample}_reverse_translated.fasta"
    log:
        "logs/{sample}_reverse_translate_proteins.log"
    conda:
        "envs/emboss.yaml"
    shell:
        """
        backtranseq -sequence {input.translated_proteins} -outseq {output.reverse_translated_fasta} 2> {log}
        """

# Step 3: KMA iterative alignment on reverse-translated CARD sequences
rule kma_iterative_alignment:
    input:
        reads="pair_merged/{sample}_merged.fastq.gz",
        reference="card_translated/{sample}_reverse_translated.fasta"
    output:
        alignment="kma_output/{sample}_aligned.fasta"
    params:
        iterations=5  # Number of iterations for KMA
    log:
        "logs/{sample}_kma_iterative.log"
    conda:
        "envs/kma_env.yaml"
    shell:
        """
        for i in $(seq 1 {params.iterations}); do
            kma -i {input.reads} -o kma_output/{wildcards.sample}_iter_$i -t 4
            cp kma_output/{wildcards.sample}_iter_$i.fsa {output.alignment}
        done 2> {log}
        """


# Step 5: Convert KMA output to FASTA for PANDAseq
rule convert_kma_to_fasta:
    input:
        kma_fasta="kma_output/{sample}_aligned.fasta"
    output:
        fasta_gz="arg_related_db/{sample}_kma_to_fasta.fasta.gz"
    log:
        "logs/{sample}_kma_to_fasta.log"
    conda:
        "envs/seqtk_env.yaml"
    shell:
        """
        seqtk seq -A {input.kma_fasta} | gzip > {output.fasta_gz} 2> {log}
        """

# Step 6: PANDAseq for ARG-related DBs
rule run_pandaseq:
    input:
        forward="arg_related_db/{sample}_R1.fastq.gz",
        reverse="arg_related_db/{sample}_R2.fastq.gz"
    output:
        merged="arg_related_db/{sample}_merged.fasta.gz"
    params:
        min_length=50,
        max_length=500,
        threshold=0.9
    threads: PANDASEQ_THREADS
    resources:
        mem_mb=PANDASEQ_MEMORY
    log:
        "logs/{sample}_pandaseq.log"
    conda:
        PANDASEQ_CONDA_ENV
    shell:
        """
        pandaseq -f {input.forward} -r {input.reverse} \
            -w >(gzip > {output.merged}) \
            -l 100 \         # Minimum length: 100bp is common for short read data
            -L 250 \         # Maximum length: 250bp works well for Illumina data
            -t 0.95 \        # Quality threshold: 0.95 is a good balance
            -o 10 \          # Minimum overlap: 10bp is standard
            -N \            # Filters sequences with Ns
            -F \            # Forward primers only
            -6 33 \         # Quality score offset (33 for Illumina)
            -A simple_bayesian \ # Algorithm choice
            -T {threads} 2> {log}
        """

# Step 7: metaSPAdes assembly
rule run_metaspades:
    input:
        merged_fasta="arg_related_db/{sample}_merged.fasta.gz"
    output:
        contigs="assembly_output/{sample}_contigs.fasta"
    log:
        "logs/{sample}_metaspades.log"
    conda:
        "envs/spades_env.yaml"
    shell:
        """
        metaspades.py -s {input.merged_fasta} -o assembly_output/{wildcards.sample} -t 8 2> {log}
        """

# Step 8: ContigExtender
rule run_contig_extender:
    input:
        contigs="assembly_output/{sample}_contigs.fasta"
    output:
        extended_contigs="contig_extender_output/{sample}_extended_contigs.fasta"
    log:
        "logs/{sample}_contig_extender.log"
    conda:
        "path/to/contig_extender_env.yaml"
    shell:
        """
        # For paired-end reads
        .dist/extender_wrapper {input.contigs} \
            --m1 {wildcards.sample}_R1_paired.fastq.gz \
            --m2 {wildcards.sample}_R2_paired.fastq.gz \
            --enable-pair-constraint \
            --extend-tolerance 3.0 \
            -o {output.extended_contigs} \
        """

# Step 8.5: filter out contigs smaller than CARD ARGs 

rule get_minimum_length_from_CARD:
    output:
        min_length=MIN_LENGTH_FILE
    script:
        "minimum_length_CARD.py"

rule filter_contigs_by_length:
    input:
        contigs=f"{OUTPUT_DIR}/extended_contigs_final.fasta",
        min_length=MIN_LENGTH_FILE
    output:
        filtered_contigs=f"{OUTPUT_DIR}/filtered_contigs.fasta"
    log:
        f"{LOG_DIR}/filter_contigs_by_length.log"
    shell:
        """
        min_len=$(cat {input.min_length})
        awk '/^>/ {{header=$0; getline seq; if (length(seq) >= ' "$min_len" ') print header "\\n" seq;}}' \
        {input.contigs} > {output.filtered_contigs}
        """


# Step 9: Taxonomic profiling with SprayNPray
rule run_spraynpray:
    input:
        extended_contigs="contig_extender_output/{sample}_extended_contigs.fasta"
    output:
        classification_report="spraynpray_output/{sample}_classification_report.txt"
    log:
        "logs/{sample}_spraynpray.log"
    conda:
        "path/to/spraynpray_env.yaml"
    shell:
        """
        spraynpray -i {input.extended_contigs} -o {output.classification_report} 2> {log}
        spray-and-pray.py -a {input.extended_contigs} -out {output.classification_report} -ref SprayNPray/db/nr.faa 2> {log}
        """

# Step 9.5: Detect structural variants (SVs) in extended contigs
rule sv_calling:
    input:
        extended_contigs="contig_extender_output/{sample}_extended_contigs.fasta"
    output:
        sv_vcf="sv_output/{sample}_structural_variants.vcf"
    log:
        "logs/{sample}_sv_calling.log"
    conda:
        "path/to/sv_calling_env.yaml"
    input:
        bam="contig_extender_output/{sample}_extended_contigs.bam",
        ref="path/to/reference.fasta"
    shell:
        """
        # First configure manta
        configManta.py --bam {input.bam} --referenceFasta {input.ref} --runDir manta_analysis
        # Then run the analysis
        manta_analysis/runWorkflow.py -j {threads} 2> {log}
        mv manta_analysis/results/variants/diploidSV.vcf.gz {output.sv_vcf}
        """

# Step 10: Identify chimeric contigs based on SVs and taxonomic classification
rule identify_chimeric_contigs:
    input:
        sv_vcf="sv_output/{sample}_structural_variants.vcf",
        kraken2_report="kraken2_output/{sample}_kraken2_report.txt",
        spraynpray_report="spraynpray_output/{sample}_classification_report.txt"
    output:
        chimeric_contigs="chimeric_output/{sample}_chimeric_contigs.fasta"
    log:
        "logs/{sample}_identify_chimeric_contigs.log"
    script:
        "scripts/identify_chimeric_contigs.py"

# Step 11: Compare chimeric contigs to plasmid and phage reference databases
rule classify_chimeric_with_plasmid_phage:
    input:
        chimeric_contigs="chimeric_output/{sample}_chimeric_contigs.fasta"
    output:
        plasmid_results="plasmid_phage_output/{sample}_plasmid_hits.txt",
        phage_results="plasmid_phage_output/{sample}_phage_hits.txt",
        combined_results="plasmid_phage_output/{sample}_classification_report.txt"
    params:
        plasmid_db="path/to/plasmid_db",
        phage_db="path/to/phage_db"
    log:
        "logs/{sample}_classify_chimeric_plasmid_phage.log"
    conda:
        "envs/blast_env.yaml"
    shell:
        """
        # Search against plasmid database
        blastn -db {params.plasmid_db} -query {input.chimeric_contigs} -out {output.plasmid_results} -num_threads 4 2>> {log}
        
        # Search against phage database
        blastn -db {params.phage_db} -query {input.chimeric_contigs} -out {output.phage_results} -num_threads 4 2>> {log}
        
        # Combine results
        cat {output.plasmid_results} {output.phage_results} > {output.combined_results}
        """

# Reverse Pipeline Rules

rule analyze_ARG_content:
    input:
        filtered_contigs=f"{OUTPUT_DIR}/filtered_contigs.fasta"
    output:
        rgi_output_dir=RGI_OUTPUT_DIR
    params:
        rgi_version="5.1.0",
        card_version="3.0.6",
        rgi_database=f"{CARD_DB_DIR}/CARD_db"
    log:
        f"{LOG_DIR}/rgi_analysis.log"
    resources:
        threads=RGI_THREADS,
        mem_mb=RGI_MEMORY
    conda:
        RGI_CONDA_ENV
    shell:
        """
        rgi main \
            --input_sequence {input.filtered_contigs} \
            --output_file {output.rgi_output_dir}/rgi_results \
            --input_type contig \
            --data {params.rgi_database} \
            --local \
            --clean \
            --threads {resources.threads} \
            --version {params.rgi_version} \
            --card_annotation_version {params.card_version} \
            2> {log}
        """

rule filter_ARG_contigs:
    input:
        rgi_output=f"{RGI_OUTPUT_DIR}/rgi_results.txt",
        filtered_contigs=f"{OUTPUT_DIR}/filtered_contigs.fasta"
    output:
        arg_contigs=f"{OUTPUT_DIR}/arg_contigs.fasta"
    params:
        cutoff=["perfect", "strict"]
    log:
        f"{LOG_DIR}/filter_ARG_contigs.log"
    script:
        "scripts/filter_ARG_contigs.py"

rule reclassify_ARG_contigs:
    input:
        arg_contigs=f"{OUTPUT_DIR}/arg_contigs.fasta",
        kraken2_db="Kraken/k2_db"
    output:
        kraken2_report=f"{LOG_DIR}/final_kraken2_report.txt",
        kraken2_output=f"{OUTPUT_DIR}/final_kraken2_output.tsv",
        classified_contigs=f"{OUTPUT_DIR}/classified_contigs.fasta",
        unclassified_contigs=f"{OUTPUT_DIR}/unclassified_contigs.fasta",
        high_confidence_contigs=f"{OUTPUT_DIR}/high_confidence_contigs.fasta",
        low_confidence_contigs=f"{OUTPUT_DIR}/low_confidence_contigs.fasta"
    params:
        kmer_length=35,
        confidence_threshold=0.8
    log:
        f"{LOG_DIR}/arg_contigs_kraken2_classification.log"
    resources:
        threads=KRAKEN2_THREADS,
        mem_mb=KRAKEN2_MEMORY
    conda:
        KRAKEN2_CONDA_ENV
    shell:
        """
        # Run Kraken2 classification
        kraken2 --db {input.kraken2_db} \
                --threads {resources.threads} \
                --output {output.kraken2_output} \
                --report {output.kraken2_report} \
                --minimum-hit-groups 2 \
                --kmer-len {params.kmer_length} \
                --classified-out {output.classified_contigs} \
                --unclassified-out {output.unclassified_contigs} \
                {input.arg_contigs} \
                2> {log}
                
        # Split contigs based on confidence score
        awk -v threshold={params.confidence_threshold} '
            /^>/ {{inheader=1; header=$0}} 
            !inheader {{seq=$0; 
                if($2 >= threshold) {{
                    print header "\\n" seq > "{output.high_confidence_contigs}"
                }} else {{
                    print header "\\n" seq > "{output.low_confidence_contigs}"
                }}
            }}' {output.kraken2_output}
        """


# Step: Evaluate contig quality (N50, L50, GC Content, Coverage)
rule evaluate_contig_quality:
    input:
        arg_contigs=f"{OUTPUT_DIR}/arg_contigs.fasta",  # The contigs FASTA file
        bam_file=f"{OUTPUT_DIR}/mapped_reads.bam"       # BAM file for calculating coverage
    output:
        quality_report=f"{OUTPUT_DIR}/contig_quality_report.txt",  # The main quality report
        length_distribution=f"{OUTPUT_DIR}/contig_length_distribution.txt"  # Contig length distribution
    log:
        f"{LOG_DIR}/evaluate_contig_quality.log"
    conda:
        "path/to/quality_metrics_env.yaml"  # Conda environment that includes BioPython and pysam
    script:
        "scripts/calculate_contig_quality.py"  # Your updated Python script
    


# Step: Map reads to all contigs
rule map_reads_to_all_contigs:
    input:
        contigs=f"{OUTPUT_DIR}/all_contigs.fasta",
        reads_r1=f"{READS_DIR}/{sample}_R1.fastq.gz",
        reads_r2=f"{READS_DIR}/{sample}_R2.fastq.gz"
    output:
        bam=f"{OUTPUT_DIR}/{sample}_mapped_reads.bam"
    log:
        f"{LOG_DIR}/map_reads_to_all_contigs.log"
    conda:
        "path/to/bwa_env.yaml"
    shell:
        """
        bwa mem -t 8 {input.contigs} {input.reads_r1} {input.reads_r2} | \
        samtools view -Sb - | \
        samtools sort -o {output.bam} 2> {log}
        """

# Step: Calculate coverage across ARG-contigs and other contigs
rule calculate_coverage_comparison:
    input:
        bam=f"{OUTPUT_DIR}/{sample}_mapped_reads.bam",
        arg_contigs=f"{OUTPUT_DIR}/arg_contigs.fasta"
    output:
        coverage_report=f"{OUTPUT_DIR}/{sample}_coverage_comparison.txt"
    log:
        f"{LOG_DIR}/calculate_coverage_comparison.log"
    conda:
        "path/to/samtools_env.yaml"
    shell:
        """
        # Calculate coverage for ARG-contigs
        samtools depth -a -b {input.arg_contigs} {input.bam} > arg_contig_coverage.txt
        
        # Calculate coverage for all contigs
        samtools depth -a {input.bam} > all_contig_coverage.txt

        # Compare the coverage between ARG-contigs and other contigs (can use a custom script)
        python scripts/compare_coverage.py arg_contig_coverage.txt all_contig_coverage.txt > {output.coverage_report}
        """

# Step: Evaluate mapping quality (for misassemblies)
rule evaluate_mapping_quality:
    input:
        bam=f"{OUTPUT_DIR}/{sample}_mapped_reads.bam"
    output:
        flagstat=f"{OUTPUT_DIR}/{sample}_flagstat.txt",
        alignment_stats=f"{OUTPUT_DIR}/{sample}_alignment_stats.txt"
    log:
        f"{LOG_DIR}/evaluate_mapping_quality.log"
    conda:
        "path/to/samtools_env.yaml"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        samtools stats {input.bam} > {output.alignment_stats}
        """

# Step: Detect coverage changes along contigs
rule detect_coverage_changes:
    input:
        bam=f"{OUTPUT_DIR}/{sample}_mapped_reads.bam",
        contigs=f"{OUTPUT_DIR}/chimeric_contigs.fasta"
    output:
        coverage=f"{OUTPUT_DIR}/{sample}_contig_coverage.txt",
        coverage_plot=f"{OUTPUT_DIR}/{sample}_coverage_plot.png"
    log:
        f"{LOG_DIR}/detect_coverage_changes.log"
    conda:
        "path/to/samtools_env.yaml"
    shell:
        """
        # Calculate coverage across contigs
        samtools depth {input.bam} > {output.coverage} 2> {log}
        
        # Custom script to generate coverage plot and check for misassemblies
        python scripts/plot_coverage.py {output.coverage} {output.coverage_plot}
        """

rule bwa_mem_ungapped:
    input:
        ref="path/to/reference.fasta",
        reads_r1="path/to/reads_R1.fastq.gz",
        reads_r2="path/to/reads_R2.fastq.gz"
    output:
        sam="path/to/alignment_k21_ungapped.sam"
    log:
        "logs/bwa_mem_ungapped.log"
    conda:
        "envs/bwa.yaml"
    shell:
        """
       bwa mem -k 21 -B 5 -O 1000 -E 1000 {input.ref} {input.reads_r1} {input.reads_r2} > {output.sam}21 2> {log} 
       bwa mem -k 31 -B 4 -O 1000 -E 1000 {input.ref} {input.reads_r1} {input.reads_r2} > {output.sam}31 2> {log}
       bwa mem -k 51 -B 3 -O 1000 -E 1000 {input.ref} {input.reads_r1} {input.reads_r2} > {output.sam}51 2> {log}
        """

rule bowtie2_ungapped:
    input:
        ref="path/to/reference",
        reads_r1="path/to/reads_R1.fastq.gz",
        reads_r2="path/to/reads_R2.fastq.gz"
    output:
        sam="path/to/alignment_bowtie2_ungapped.sam"
    log:
        "logs/bowtie2_ungapped.log"
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -x {input.ref} -1 {input.reads_r1} -2 {input.reads_r2} \
            -S {output.sam} --score-min L,0,-0.05 --mp 4,2 --np 1 --rdg 1000,1000 \
            --rfg 1000,1000 --no-1mm-upfront --very-sensitive 2> {log}
        """


rule kma_ungapped:
    input:
        reads="path/to/reads.fastq",
        ref_db="path/to/kma_database"
    output:
        output_prefix="path/to/kma_ungapped_alignment"
    log:
        "logs/kma_ungapped.log"
    conda:
        "envs/kma.yaml"
    shell:
        """
        kma -i {input.reads} -o {output.output_prefix} -t_db {input.ref_db} \
            -ID 0.95 -mem_mode 1 -ef 1 -g 0 2> {log}
        """

rule minimap2_ungapped:
    input:
        ref="path/to/reference.fasta",
        reads="path/to/reads.fastq"
    output:
        sam="path/to/alignment_minimap2_ungapped.sam"
    log:
        "logs/minimap2_ungapped.log"
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        minimap2 -x sr -a {input.ref} {input.reads} -A 2 -B 5 -O 1000,1000 -E 1000,1000 > {output.sam} 2> {log}
        """

rule tblastn_mapping:
    input:
        protein_seqs="path/to/protein_sequences.faa",
        ref_db="path/to/reference_nucleotide.fasta"
    output:
        tblastn_output="path/to/tblastn_output.tsv"
    log:
        "logs/tblastn_mapping.log"
    conda:
        "envs/blast.yaml"
    shell:
        """
        tblastn -query {input.protein_seqs} -db {input.ref_db} \
            -outfmt 6 -evalue 1e-5 -out {output.tblastn_output} \
            -gapopen 5 -gapextend 2 -matrix BLOSUM62 2> {log}
        """

rule blastp_mapping:
    input:
        protein_seqs="path/to/protein_sequences.faa",
        ref_db="path/to/reference_protein.faa"
    output:
        blastp_output="path/to/blastp_output.tsv"
    log:
        "logs/blastp_mapping.log"
    conda:
        "envs/blast.yaml"
    shell:
        """
        blastp -query {input.protein_seqs} -db {input.ref_db} \
            -outfmt 6 -evalue 1e-5 -out {output.blastp_output} \
            -gapopen 5 -gapextend 2 -matrix BLOSUM62 2> {log}
        """


rule mash_sketch:
    input:
        fasta="path/to/contigs_or_plasmids.fasta"
    output:
        sketch="path/to/sketch.msh"
    shell:
        """
        mash sketch -o {output.sketch} {input.fasta}
        """

rule bray_curtis_diversity:
    input:
        abundance_matrix="path/to/abundance_table.csv"
    output:
        bray_curtis_matrix="path/to/bray_curtis_distances.csv",
        pcoa_plot="path/to/pcoa_plot.png"
    script:
        "scripts/calculate_bray_curtis.py"  # Or an R script


rule mash_sketch_mags:
    input:
        fasta="path/to/MAGs.fasta"
    output:
        sketch="path/to/mag_sketch.msh"
    shell:
        """
        mash sketch -o {output.sketch} {input.fasta}
        """

rule bray_curtis_diversity_mags:
    input:
        presence_absence_matrix="path/to/presence_absence_matrix.csv"
    output:
        bray_curtis_matrix="path/to/bray_curtis_distances.csv",
        pcoa_plot="path/to/pcoa_plot.png"
    script:
        "scripts/calculate_bray_curtis_presence_absence.py"  # Or an R script


##Phages?
rule extract_phage_genes:
    input:
        chimeric_contigs="path/to/chimeric_contigs.fasta"
    output:
        phage_genes="path/to/phage_genes.fasta"
    shell:
        """
        # Your custom script to extract phage-associated genes from the chimeric contigs
        extract_phage_genes.py {input.chimeric_contigs} > {output.phage_genes}
        """
rule build_phylogenetic_tree:
    input:
        phage_genes="path/to/phage_genes.fasta"
    output:
        tree="path/to/phage_tree.newick"
    log:
        "logs/phylogenetic_tree.log"
    conda:
        "envs/phylogenetics.yaml"
    shell:
        """
        fasttree -nt {input.phage_genes} > {output.tree} 2> {log}
        """
rule visualize_phylogenetic_tree:
    input:
        tree="path/to/phage_tree.newick"
    output:
        visualization="path/to/tree_visualization.png"
    shell:
        """
        # Example using iTOL API or FigTree for visualization
        visualize_tree.py {input.tree} > {output.visualization}
        """


# PlasPredict Pipeline

rule run_plaspredict:
    input:
        contigs=f"{OUTPUT_DIR}/arg_contigs.fasta"
    output:
        plaspredict_output_dir=PLASPREDICT_OUTPUT_DIR
    params:
        plaspredict_env="plaspredict_env"
    log:
        f"{LOG_DIR}/run_plaspredict.log"
    resources:
        threads=PLASPREDICT_THREADS,
        mem_mb=PLASPREDICT_MEMORY
    conda:
        PLASPREDICT_CONDA_ENV
    shell:
        """
        plaspredict --input {input.contigs} --output {output.plaspredict_output_dir} \
        --threads {resources.threads} 2> {log}
        """

# MGE Identification and Categorization

rule hmm_scan_MGEs:
    input:
        predicted_proteins=PROTEIN_OUTPUT,
        pfam_db=PFAM_DB,
        tnppred_db=TNPPRED_DB
    output:
        pfam_hmm_output=f"{OUTPUT_DIR}/pfam_hmm_output.txt",
        tnppred_hmm_output=f"{OUTPUT_DIR}/tnppred_hmm_output.txt"
    log:
        f"{LOG_DIR}/hmm_scan_MGEs.log"
    resources:
        threads=HMMER_THREADS,
        mem_mb=HMMER_MEMORY
    conda:
        HMMER_CONDA_ENV
    shell:
        """
        hmmsearch --cpu {resources.threads} \
                  --domtblout {output.pfam_hmm_output} \
                  {input.pfam_db} \
                  {input.predicted_proteins} \
                  2>> {log}

        hmmsearch --cpu {resources.threads} \
                  --domtblout {output.tnppred_hmm_output} \
                  {input.tnppred_db} \
                  {input.predicted_proteins} \
                  2>> {log}
        """

rule filter_and_assign_hits:
    input:
        pfam_hmm_output=f"{OUTPUT_DIR}/pfam_hmm_output.txt",
        tnppred_hmm_output=f"{OUTPUT_DIR}/tnppred_hmm_output.txt"
    output:
        best_hits=BEST_HITS
    params:
        e_value_threshold=1e-5,
        orf_distance=10
    log:
        f"{LOG_DIR}/filter_and_assign_hits.log"
    script:
        "scripts/filter_and_assign_hits.py"

rule categorize_MGE_domains:
    input:
        best_hits=BEST_HITS,
        rgi_output=f"{RGI_OUTPUT_DIR}/rgi_results.txt"
    output:
        categorized_MGEs=f"{OUTPUT_DIR}/categorized_MGEs.txt"
    log:
        f"{LOG_DIR}/categorize_MGE_domains.log"
    script:
        "scripts/categorize_MGE_domains.py"

