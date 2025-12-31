import os
import random

# Store generated parameters in a dictionary to avoid repetition
generated_parameters = {}

import random

# Store generated parameters in a dictionary to avoid repetition
generated_parameters = {}

def random_trimmomatic_parameters(log_file, iteration, sample):
    # Use the iteration as the key for seed
    seed_key = (iteration, 'trimmomatic_seed')
    if seed_key not in generated_parameters:
        # Generate a single seed per iteration
        seed = random.randint(0, 100)
        generated_parameters[seed_key] = seed
    
    # Use the combination of iteration and sample for other parameters
    param_key = (iteration, sample, 'trimmomatic')
    if param_key not in generated_parameters:
        minlen = random.randint(30, 50)
        swindow = random.randint(4, 10)
        squality = random.randint(15, 30)
        generated_parameters[param_key] = (generated_parameters[seed_key], minlen, swindow, squality)
        
        # Save the parameters to a TSV file only once
        with open(log_file, 'a') as f:
            f.write(f"Trimmomatic\t{iteration}\t{sample}\t{generated_parameters[seed_key]}\t{minlen}\t{swindow}\t{squality}\n")
    return generated_parameters[param_key]

def random_fastp_parameters(log_file, iteration, sample):
    # Use the iteration as the key for seed (if needed)
    seed_key = (iteration, 'fastp_seed')
    if seed_key not in generated_parameters:
        # Generate a single seed per iteration (if needed)
        seed = random.randint(0, 100)
        generated_parameters[seed_key] = seed
    
    param_key = (iteration, sample, 'fastp')
    if param_key not in generated_parameters:
        quality_cut = random.randint(20, 30)
        length_cut = random.randint(30, 50)
        n_base_limit = random.randint(0, 10)
        generated_parameters[param_key] = (generated_parameters[seed_key], quality_cut, length_cut, n_base_limit)
        
        # Save the parameters to a TSV file only once
        with open(log_file, 'a') as f:
            f.write(f"Fastp\t{iteration}\t{sample}\t{generated_parameters[seed_key]}\t{quality_cut}\t{length_cut}\t{n_base_limit}\n")
    return generated_parameters[param_key]


def random_cutadapt_parameters(log_file, iteration, sample):
    # Use the iteration as the key for seed (if needed)
    seed_key = (iteration, 'cutadapt_seed')
    if seed_key not in generated_parameters:
        # Generate a single seed per iteration (if needed)
        seed = random.randint(0, 100)
        generated_parameters[seed_key] = seed
    
    param_key = (iteration, sample, 'cutadapt')
    if param_key not in generated_parameters:
        adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  # Adapter sequence for Read 1
        adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  # Adapter sequence for Read 2
        error_rate = round(random.uniform(0.1, 0.2), 2)
        minimum_length = random.randint(25, 40)
        overlap = random.randint(3, 10)
        generated_parameters[param_key] = (generated_parameters[seed_key], adapter1, adapter2, error_rate, minimum_length, overlap)
        
        # Save the parameters to a TSV file only once
        with open(log_file, 'a') as f:
            f.write(f"Cutadapt\t{iteration}\t{sample}\t{generated_parameters[seed_key]}\t{adapter1}\t{adapter2}\t{error_rate}\t{minimum_length}\t{overlap}\n")
    return generated_parameters[param_key]


def random_bbduk_parameters(log_file, iteration, sample):
    # Use the iteration as the key for seed (if needed)
    seed_key = (iteration, 'bbduk_seed')
    if seed_key not in generated_parameters:
        # Generate a single seed per iteration (if needed)
        seed = random.randint(0, 100)
        generated_parameters[seed_key] = seed
    
    param_key = (iteration, sample, 'bbduk')
    if param_key not in generated_parameters:
        ktrim = random.choice(['r', 'l'])  # 'r' for right end, 'l' for left end
        k = random.randint(19, 27)  # k-mer length
        mink = random.randint(8, 12)  # minimum k-mer length for right end trimming
        hdist = random.randint(0, 2)  # hamming distance for k-mer matching
        minlen = random.randint(30, 50)  # minimum length of reads after trimming
        qtrim = random.choice(['rl', 'f', 'r', 'l'])  # quality trimming mode
        trimq = random.randint(10, 20)  # quality threshold for trimming
        minkmerhits = random.randint(1, 2)  # Minimum k-mer hits required to consider a read as matching
        minkmerfraction = round(random.uniform(0.0, 0.5), 2)  # Minimum k-mer fraction
        generated_parameters[param_key] = (generated_parameters[seed_key], ktrim, k, mink, hdist, minlen, qtrim, trimq, minkmerhits, minkmerfraction)
        
        # Save the parameters to a TSV file only once
        with open(log_file, 'a') as f:
            f.write(f"BBDuk\t{iteration}\t{sample}\t{generated_parameters[seed_key]}\t{ktrim}\t{k}\t{mink}\t{hdist}\t{minlen}\t{qtrim}\t{trimq}\t{minkmerhits}\t{minkmerfraction}\n")
    return generated_parameters[param_key]


def random_sickle_parameters(log_file, iteration, sample):
    # Use the iteration as the key for seed (if needed)
    seed_key = (iteration, 'sickle_seed')
    if seed_key not in generated_parameters:
        # Generate a single seed per iteration (if needed)
        seed = random.randint(0, 100)
        generated_parameters[seed_key] = seed
    
    param_key = (iteration, sample, 'sickle')
    if param_key not in generated_parameters:
        quality_threshold = random.randint(20, 30)  # quality threshold
        length_threshold = random.randint(25, 50)  # length threshold after trimming
        generated_parameters[param_key] = (generated_parameters[seed_key], quality_threshold, length_threshold)
        
        # Save the parameters to a TSV file only once
        with open(log_file, 'a') as f:
            f.write(f"Sickle\t{iteration}\t{sample}\t{generated_parameters[seed_key]}\t{quality_threshold}\t{length_threshold}\n")
    return generated_parameters[param_key]


ITERATIONS = 3  # Number of times to loop for each file


rule all:
    input:
        expand("output_dir/fastp_output/iteration_{iteration}/{sample}_R1_fastp_trimmed.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R1', '') for file in os.listdir("raw_reads") if file.endswith('_R1.fastq.gz')]),
        expand("output_dir/fastp_output/iteration_{iteration}/{sample}_R2_fastp_trimmed.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R2', '') for file in os.listdir("raw_reads") if file.endswith('_R2.fastq.gz')]),
        expand("output_dir/trimmomatic_output/iteration_{iteration}/{sample}_R1_paired.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R1', '') for file in os.listdir("raw_reads") if file.endswith('_R1.fastq.gz')]),
        expand("output_dir/trimmomatic_output/iteration_{iteration}/{sample}_R1_unpaired.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R2', '') for file in os.listdir("raw_reads") if file.endswith('_R2.fastq.gz')]),
        expand("output_dir/trimmomatic_output/iteration_{iteration}/{sample}_R2_paired.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R1', '') for file in os.listdir("raw_reads") if file.endswith('_R1.fastq.gz')]),
        expand("output_dir/trimmomatic_output/iteration_{iteration}/{sample}_R2_unpaired.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R2', '') for file in os.listdir("raw_reads") if file.endswith('_R2.fastq.gz')]),
        expand("output_dir/cutadapt_output/iteration_{iteration}/{sample}_R1_cutadapt_trimmed.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R1', '') for file in os.listdir("raw_reads") if file.endswith('_R1.fastq.gz')]),
        expand("output_dir/cutadapt_output/iteration_{iteration}/{sample}_R2_cutadapt_trimmed.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R2', '') for file in os.listdir("raw_reads") if file.endswith('_R2.fastq.gz')]),
        expand("output_dir/bbduk_output/iteration_{iteration}/{sample}_R1_bbduk_trimmed.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R1', '') for file in os.listdir("raw_reads") if file.endswith('_R1.fastq.gz')]),
        expand("output_dir/bbduk_output/iteration_{iteration}/{sample}_R2_bbduk_trimmed.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R2', '') for file in os.listdir("raw_reads") if file.endswith('_R2.fastq.gz')]),
        expand("output_dir/sickle_output/iteration_{iteration}/{sample}_R1_sickle_trimmed.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R1', '') for file in os.listdir("raw_reads") if file.endswith('_R1.fastq.gz')]),
        expand("output_dir/sickle_output/iteration_{iteration}/{sample}_R2_sickle_trimmed.fastq.gz", 
               iteration=range(1, ITERATIONS + 1), 
               sample=[os.path.splitext(os.path.splitext(file)[0])[0].replace('_R2', '') for file in os.listdir("raw_reads") if file.endswith('_R2.fastq.gz')])




### trimmomatic
LOG_FILE_trimmomatic = "trimmomatic_params.tsv"

# Ensure the log file has headers only once
if not os.path.exists(LOG_FILE_trimmomatic):
    with open(LOG_FILE_trimmomatic, 'w') as f:
        f.write("Tool\tIteration\tSample\tSeed\tMinLen\tSlidingWindow\n")


rule trimmomatic:
    input:
        forward_reads="raw_reads/{sample}_R1.fastq.gz",
        reverse_reads="raw_reads/{sample}_R2.fastq.gz"
    output:
        forward_paired="output_dir/trimmomatic_output/iteration_{iteration}/{sample}_R1_paired.fastq.gz",
        forward_unpaired="output_dir/trimmomatic_output/iteration_{iteration}/{sample}_R1_unpaired.fastq.gz",
        reverse_paired="output_dir/trimmomatic_output/iteration_{iteration}/{sample}_R2_paired.fastq.gz",
        reverse_unpaired="output_dir/trimmomatic_output/iteration_{iteration}/{sample}_R2_unpaired.fastq.gz"
    params:
        log_file=LOG_FILE_trimmomatic,
        seed=lambda wildcards: random_trimmomatic_parameters(LOG_FILE_trimmomatic, wildcards.iteration, wildcards.sample)[0],
        minlen=lambda wildcards: random_trimmomatic_parameters(LOG_FILE_trimmomatic, wildcards.iteration, wildcards.sample)[1],
        swindow=lambda wildcards: random_trimmomatic_parameters(LOG_FILE_trimmomatic, wildcards.iteration, wildcards.sample)[2],
        squality=lambda wildcards: random_trimmomatic_parameters(LOG_FILE_trimmomatic, wildcards.iteration, wildcards.sample)[3]
    shell:
        """
        for i in $(seq 0 $(expr $(echo {input.forward_reads} | wc -w) - 1)); do \
            trimmomatic PE -phred33 \
                $(echo {input.forward_reads} | cut -d ' ' -f $((i+1))) \
                $(echo {input.reverse_reads} | cut -d ' ' -f $((i+1))) \
                $(echo {output.forward_paired} | cut -d ' ' -f $((i+1))) \
                $(echo {output.forward_unpaired} | cut -d ' ' -f $((i+1))) \
                $(echo {output.reverse_paired} | cut -d ' ' -f $((i+1))) \
                $(echo {output.reverse_unpaired} | cut -d ' ' -f $((i+1))) \
                ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:1:true \
                SLIDINGWINDOW:{params.swindow}:{params.squality} MINLEN:{params.minlen}; \
        done
        """



### fastp
LOG_FILE_fastp = "fastp_params.tsv"

# Ensure the log file has headers only once
if not os.path.exists(LOG_FILE_fastp):
    with open(LOG_FILE_fastp, 'w') as f:
        f.write("Tool\tIteration\tSample\tSeed\tQualityCut\tLengthCut\tNBaseLimit\n")

        
rule fastp:
    input:
        forward_paired="raw_reads/{sample}_R1.fastq.gz",
        reverse_paired="raw_reads/{sample}_R2.fastq.gz"
    output:
        forward_trimmed="output_dir/fastp_output/iteration_{iteration}/{sample}_R1_fastp_trimmed.fastq.gz",
        reverse_trimmed="output_dir/fastp_output/iteration_{iteration}/{sample}_R2_fastp_trimmed.fastq.gz"
    params:
        log_file=LOG_FILE_fastp,
        quality_cut=lambda wildcards: random_fastp_parameters(LOG_FILE_fastp, wildcards.iteration, wildcards.sample)[0],
        length_cut=lambda wildcards: random_fastp_parameters(LOG_FILE_fastp, wildcards.iteration, wildcards.sample)[1],
        n_base_limit=lambda wildcards: random_fastp_parameters(LOG_FILE_fastp, wildcards.iteration, wildcards.sample)[2]
    shell:
        """
        for i in $(seq 0 $(expr $(echo {input.forward_paired} | wc -w) - 1)); do \
            fastp -i $(echo {input.forward_paired} | cut -d ' ' -f $((i+1))) \
                  -I $(echo {input.reverse_paired} | cut -d ' ' -f $((i+1))) \
                  -o $(echo {output.forward_trimmed} | cut -d ' ' -f $((i+1))) \
                  -O $(echo {output.reverse_trimmed} | cut -d ' ' -f $((i+1))) \
                  --qualified_quality_phred {params.quality_cut} \
                  --length_required {params.length_cut} \
                  --n_base_limit {params.n_base_limit}; \
        done
        """



### Cutadapt
LOG_FILE_cutadapt = "cutadapt_params.tsv"

# Ensure the log file has headers only once
if not os.path.exists(LOG_FILE_cutadapt):
    with open(LOG_FILE_cutadapt, 'w') as f:
        f.write("Tool\tIteration\tSample\tSeed\tAdapter1\tAdapter2\tErrorRate\tMinimumLength\tOverlap\n")

rule cutadapt:
    input:
        forward_paired="raw_reads/{sample}_R1.fastq.gz",
        reverse_paired="raw_reads/{sample}_R2.fastq.gz"
    output:
        forward_trimmed="output_dir/cutadapt_output/iteration_{iteration}/{sample}_R1_cutadapt_trimmed.fastq.gz",
        reverse_trimmed="output_dir/cutadapt_output/iteration_{iteration}/{sample}_R2_cutadapt_trimmed.fastq.gz"
    params:
        log_file=LOG_FILE_cutadapt,
        adapter1=lambda wildcards: random_cutadapt_parameters(LOG_FILE_cutadapt, wildcards.iteration, wildcards.sample)[1],
        adapter2=lambda wildcards: random_cutadapt_parameters(LOG_FILE_cutadapt, wildcards.iteration, wildcards.sample)[2],
        error_rate=lambda wildcards: random_cutadapt_parameters(LOG_FILE_cutadapt, wildcards.iteration, wildcards.sample)[3],
        minimum_length=lambda wildcards: random_cutadapt_parameters(LOG_FILE_cutadapt, wildcards.iteration, wildcards.sample)[4],
        overlap=lambda wildcards: random_cutadapt_parameters(LOG_FILE_cutadapt, wildcards.iteration, wildcards.sample)[5]
    shell:
        """
        for i in $(seq 0 $(expr $(echo {input.forward_paired} | wc -w) - 1)); do \
            cutadapt -a {params.adapter1} \
                    -A {params.adapter2} \
                    -o $(echo {output.forward_trimmed} | cut -d ' ' -f $((i+1))) \
                    -p $(echo {output.reverse_trimmed} | cut -d ' ' -f $((i+1))) \
                    --error-rate={params.error_rate} \
                    --minimum-length={params.minimum_length} \
                    --overlap={params.overlap} \
                    $(echo {input.forward_paired} | cut -d ' ' -f $((i+1))) \
                    $(echo {input.reverse_paired} | cut -d ' ' -f $((i+1))); \
        done
        """


#### BBDUK
LOG_FILE_bbduk = "bbduk_params.tsv"

# Ensure the log file has headers only once
if not os.path.exists(LOG_FILE_bbduk):
    with open(LOG_FILE_bbduk, 'w') as f:
        f.write("Tool\tIteration\tSample\tSeed\tKTrim\tK\tMinK\tHDist\tMinLen\tQTrim\tTrimQ\tMinKmerHits\tMinKmerFraction\n")

rule bbduk:
    input:
        forward_paired="raw_reads/{sample}_R1.fastq.gz",
        reverse_paired="raw_reads/{sample}_R2.fastq.gz"
    output:
        forward_trimmed="output_dir/bbduk_output/iteration_{iteration}/{sample}_R1_bbduk_trimmed.fastq.gz",
        reverse_trimmed="output_dir/bbduk_output/iteration_{iteration}/{sample}_R2_bbduk_trimmed.fastq.gz"
    threads: 8
    params:
        log_file=LOG_FILE_bbduk,
        ktrim=lambda wildcards: random_bbduk_parameters(LOG_FILE_bbduk, wildcards.iteration, wildcards.sample)[1],
        k=lambda wildcards: random_bbduk_parameters(LOG_FILE_bbduk, wildcards.iteration, wildcards.sample)[2],
        mink=lambda wildcards: random_bbduk_parameters(LOG_FILE_bbduk, wildcards.iteration, wildcards.sample)[3],
        hdist=lambda wildcards: random_bbduk_parameters(LOG_FILE_bbduk, wildcards.iteration, wildcards.sample)[4],
        minlen=lambda wildcards: random_bbduk_parameters(LOG_FILE_bbduk, wildcards.iteration, wildcards.sample)[5],
        qtrim=lambda wildcards: random_bbduk_parameters(LOG_FILE_bbduk, wildcards.iteration, wildcards.sample)[6],
        trimq=lambda wildcards: random_bbduk_parameters(LOG_FILE_bbduk, wildcards.iteration, wildcards.sample)[7],
        minkmerhits=lambda wildcards: random_bbduk_parameters(LOG_FILE_bbduk, wildcards.iteration, wildcards.sample)[8],
        minkmerfraction=lambda wildcards: random_bbduk_parameters(LOG_FILE_bbduk, wildcards.iteration, wildcards.sample)[9],
        ref="adapters"  # Assuming adapters are the common contaminants you're removing
    shell:
        """
        for i in $(seq 0 $(expr $(echo {input.forward_paired} | wc -w) - 1)); do \
            bbduk.sh in1=$(echo {input.forward_paired} | cut -d ' ' -f $((i+1))) \
                     in2=$(echo {input.reverse_paired} | cut -d ' ' -f $((i+1))) \
                     out1=$(echo {output.forward_trimmed} | cut -d ' ' -f $((i+1))) \
                     out2=$(echo {output.reverse_trimmed} | cut -d ' ' -f $((i+1))) \
                     ktrim={params.ktrim} k={params.k} mink={params.mink} \
                     hdist={params.hdist} minlen={params.minlen} \
                     qtrim={params.qtrim} trimq={params.trimq} \
                     minkmerhits={params.minkmerhits} minkmerfraction={params.minkmerfraction} \
                     maskMiddle=f \
                     literal=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -da; \
        done
        """



LOG_FILE_sickle = "sickle_params.tsv"

# Ensure the log file has headers only once
if not os.path.exists(LOG_FILE_sickle):
    with open(LOG_FILE_sickle, 'w') as f:
        f.write("Tool\tIteration\tSample\tSeed\tQualityThreshold\tLengthThreshold\n")

rule sickle:
    input:
        forward_paired="raw_reads/{sample}_R1.fastq.gz",
        reverse_paired="raw_reads/{sample}_R2.fastq.gz"
    output:
        forward_trimmed="output_dir/sickle_output/iteration_{iteration}/{sample}_R1_sickle_trimmed.fastq.gz",
        reverse_trimmed="output_dir/sickle_output/iteration_{iteration}/{sample}_R2_sickle_trimmed.fastq.gz",
        single_trimmed="output_dir/sickle_output/iteration_{iteration}/{sample}_singles_sickle_trimmed.fastq.gz"
    params:
        log_file=LOG_FILE_sickle,
        quality_threshold=lambda wildcards: random_sickle_parameters(LOG_FILE_sickle, wildcards.iteration, wildcards.sample)[1],
        length_threshold=lambda wildcards: random_sickle_parameters(LOG_FILE_sickle, wildcards.iteration, wildcards.sample)[2]
    shell:
        """
        for i in $(seq 0 $(expr $(echo {input.forward_paired} | wc -w) - 1)); do \
            sickle pe -f $(echo {input.forward_paired} | cut -d ' ' -f $((i+1))) \
                      -r $(echo {input.reverse_paired} | cut -d ' ' -f $((i+1))) \
                      -t sanger \
                      -o $(echo {output.forward_trimmed} | cut -d ' ' -f $((i+1))) \
                      -p $(echo {output.reverse_trimmed} | cut -d ' ' -f $((i+1))) \
                      -s $(echo {output.single_trimmed} | cut -d ' ' -f $((i+1))) \
                      -q {params.quality_threshold} \
                      -l {params.length_threshold} \
                      -g; \
        done
        """
       
