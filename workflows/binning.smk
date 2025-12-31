# Environment settings for binning tools
METABAT_ENV = "metabat2_env"
MAXBIN_ENV = "maxbin2_env"
CONCOCT_ENV = "concoct_env"
DASTOOL_ENV = "dastool_env"
CHECKM2_ENV = "checkm2_env"

rule all:
    input:
        # Binning outputs moved to separate file
        expand("binning_output/{sample}_binning_output/{binner}_bins",
               sample=SAMPLES,
               binner=['metabat2', 'maxbin2', 'concoct']),
        # Quality assessment
        expand("checkm2_output/{sample}/quality_report.tsv",
               sample=SAMPLES)


rule metawrap_binning:
    input:
        contigs="scaffolds/{site}/contigextender_scaffolds.fasta",  # Updated to match contigextender output
        read1=f"{FASTUNIQ_DIR}/{{site}}_1.fastq",
        read2=f"{FASTUNIQ_DIR}/{{site}}_2.fastq"
    output:
        concoct_bins=f"{BINNING_DIR}/{{site}}_binning_output/concoct_bins.scaffolds2bin.tsv",
        metabat_bins=f"{BINNING_DIR}/{{site}}_binning_output/metabat2_bins.scaffolds2bin.tsv",
        maxbin_bins=f"{BINNING_DIR}/{{site}}_binning_output/maxbin2_bins.scaffolds2bin.tsv"
    params:
        min_contig_len=1500,
        memory=32,
        outdir=f"{BINNING_DIR}/{{site}}_binning_output"
    threads: 16
    log:
        f"{LOG_DIR}/metawrap_binning_{{site}}.log"
    conda: METAWRAP_BINNING_ENV
    shell:
        """
        mkdir -p {params.outdir}

        metawrap binning \
            -o {params.outdir} \
            -t {threads} \
            -a {input.contigs} \
            --metabat2 --maxbin2 --concoct \
            --universal \
            -m {params.memory} \
            -l {params.min_contig_len} \
            {input.read1} {input.read2} 2>&1 | tee {log}

        # Verify bin creation
        if ! ls {params.outdir}/*bins.scaffolds2bin.tsv 1> /dev/null 2>&1; then
            echo "Error: No bins were created" >> {log}
            exit 1
        fi
        """

rule run_dastool:
    input:
        scaffolds="scaffolds/{site}/contigextender_scaffolds.fasta",
        metabat_bins=f"{BINNING_DIR}/{{site}}_binning_output/metabat2_bins.scaffolds2bin.tsv",
        maxbin_bins=f"{BINNING_DIR}/{{site}}_binning_output/maxbin2_bins.scaffolds2bin.tsv",
        concoct_bins=f"{BINNING_DIR}/{{site}}_binning_output/concoct_bins.scaffolds2bin.tsv"
    output:
        dastool_out=directory(f"{DASTOOL_DIR}/DAS_Tool_bins_{{site}}")
    params:
        labels="maxbin2,metabat2,concoct",
        search_engine="diamond",
        score_threshold=0.5
    threads: 8
    conda: "dastool_env"
    log:
        "logs/dastool_{site}.log"
    shell:
        """
        DAS_Tool \
            -i {input.maxbin_bins},{input.metabat_bins},{input.concoct_bins} \
            -l {params.labels} \
            -c {input.scaffolds} \
            -o {output.dastool_out} \
            --search_engine {params.search_engine} \
            --score_threshold {params.score_threshold} \
            --write_bins 1 \
            --threads {threads} \
            2>&1 | tee {log}
        """

# Rule to reassemble bins with MetaWRAP by site
rule reassemble_bins:
    input:
        read1_clean=f"{FASTUNIQ_DIR}/{{site}}_R1_uniq.fastq",
        read2_clean=f"{FASTUNIQ_DIR}/{{site}}_R2_uniq.fastq",
        dastool_bins=f"{DASTOOL_DIR}/DAS_Tool_bins_{{site}}"
    output:
        reassembly_out=directory(f"{REASSEMBLY_DIR}/{{site}}_reassembly_output/assembly.fasta")  # Explicit output
    threads: METAWRAP_REASSEMBLY_THREADS
    log:
        f"{LOG_DIR}/reassemble_bins_{{site}}.log"
    conda: METAWRAP_REASSEMBLY_ENV
    shell:
        """
        metawrap reassemble_bins \
            -o {REASSEMBLY_DIR}/{{site}}_reassembly_output \
            -1 {input.read1_clean} \
            -2 {input.read2_clean} \
            -b {input.dastool_bins} \
            -c {COMPLETENESS} \
            -l {MIN_CONTIG_LENGTH} \
            -t {threads} \
            -m {METAWRAP_REASSEMBLY_MEMORY} &> {log}
        """

# Rule to run MAGpurify to clean the bins by site
rule run_magpurify:
    input:
        mag=f"{REASSEMBLY_DIR}/{{site}}_reassembly_output/assembly.fasta"
    output:
        cleaned_mag=f"{MAGPURIFY_DIR}/{{site}}_cleaned.fasta"
    threads: MAGPURIFY_THREADS
    log:
        f"{LOG_DIR}/run_magpurify_{{site}}.log"
    conda: MAGPURIFY_ENV
    shell:
        """
        for module in {MAGPURIFY_MODULES}; do
            magpurify $module {input.mag} {MAGPURIFY_DIR}
        done
        magpurify clean-bin {input.mag} {MAGPURIFY_DIR} {output.cleaned_mag} &> {log}
        """



# Rule to assess assembly quality with MetaQUAST by site
rule run_metaquast:
    input:
        magpurify_bins=expand("{magpurify_out}/*.fasta", magpurify_out=f"{MAGPURIFY_DIR}/{{site}}_magpurify_output")
    output:
        metaquast_out=directory(f"{METAQUAST_DIR}/{{site}}_metaquast_output")
    threads: METAQUAST_THREADS
    log:
        f"{LOG_DIR}/run_metaquast_{{site}}.log"
    conda: METAQUAST_ENV
    shell:
        """
        metaquast.py {input.magpurify_bins} \
                     -o {output.metaquast_out} \
                     -m {MIN_CONTIG_LENGTH} \
                     -t {threads} &> {log}
        """

# Rule to dereplicate bins using dRep
rule run_drep:
    input:
        magpurify_bins=expand(f"{MAGPURIFY_DIR}/{{site}}_cleaned.fasta", site=get_sites())
    output:
        drep_out=directory(f"{DREP_DIR}/dereplicated_genomes")
    threads: DREP_THREADS
    log:
        f"{LOG_DIR}/run_drep.log"
    conda: DREP_ENV
    shell:
        """
        dRep dereplicate {output.drep_out} \
            -g {" ".join(input.magpurify_bins)} \
            -p {threads} \
            -comp {COMPLETENESS} \
            -con {CONTAMINATION} \
            --S_algorithm {DREP_S_ALGORITHM} \
            -pa {DREP_PA} \
            -sa {DREP_SA} \
            --cov_thresh {DREP_COV_THRESH} \
            --clusterAlg {DREP_CLUSTER_ALG} &> {log}
        """

# Rule to run CheckM2 for bin quality estimation
rule run_checkm2:
    input:
        bins=f"{DREP_DIR}/dereplicated_genomes/dereplicated_genomes.fasta"
    output:
        checkm2_report=f"{DREP_DIR}/checkm2_output/quality_report.tsv"
    threads: 30  # Adjust based on your system resources
    conda: CHECKM2_ENV
    shell:
        """
        checkm2 predict --threads {threads} --input {input.bins} --output-directory {DREP_DIR}/checkm2_output &> {log}
        """
