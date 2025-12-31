import os
import argparse
import subprocess
import pysam
from Bio import SeqIO
import logging
import sys
from datetime import datetime

def setup_logging(output_dir):
    """Setup logging to both file and console"""
    create_directory(output_dir)
    log_file = os.path.join(output_dir, "mobilome_analysis.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger('mobilome')

def create_directory(path):
    os.makedirs(path, exist_ok=True)
    logging.info(f"Created directory: {path}")

def run_command(cmd, description=""):
    """Run shell command with logging"""
    logging.info(f"Running {description}..." if description else f"Running command: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=True, 
                              stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                              text=True)
        if result.stdout:
            logging.debug(f"Command output: {result.stdout}")
        logging.info(f"Successfully completed {description}" if description else "Command completed successfully")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with error: {e.stderr}")
        raise

def extract_arg_regions(contigs_file, blast_results, output_bed, output_fasta):
    """Extract ARG regions ±5kb and create BED/FASTA files"""
    logging.info("Starting ARG region extraction...")
    
    regions = []
    logging.info(f"Reading BLAST results from {blast_results}")
    with open(blast_results) as f:
        for i, line in enumerate(f, 1):
            contig, start, end = line.strip().split('\t')
            start = max(0, int(start) - 5000)
            end = int(end) + 5000
            regions.append((contig, start, end))
        logging.info(f"Processed {i} BLAST result regions")
    
    # Sort and merge overlapping regions
    logging.info("Sorting and merging overlapping regions...")
    regions.sort()
    merged = []
    for region in regions:
        if not merged or merged[-1][0] != region[0] or merged[-1][2] < region[1]:
            merged.append(region)
        else:
            merged[-1] = (merged[-1][0], merged[-1][1], max(merged[-1][2], region[2]))
    logging.info(f"Merged into {len(merged)} distinct regions")
    
    # Write BED file
    logging.info(f"Writing BED file to {output_bed}")
    with open(output_bed, 'w') as f:
        for contig, start, end in merged:
            f.write(f"{contig}\t{start}\t{end}\n")
    
    # Extract sequences
    logging.info(f"Extracting sequences from {contigs_file}")
    sequences_written = 0
    with open(output_fasta, 'w') as out:
        for record in SeqIO.parse(contigs_file, "fasta"):
            for contig, start, end in merged:
                if record.id == contig:
                    sequence = record.seq[start:end]
                    out.write(f">{contig}_{start}_{end}\n{sequence}\n")
                    sequences_written += 1
    
    logging.info(f"Wrote {sequences_written} sequences to {output_fasta}")

def process_mobilome(reads1, reads2, regions_fasta, output_dir, min_coverage):
    """Process all mobilome analysis steps using KMA instead of BWA"""
    create_directory(output_dir)
    
    # Output files
    kma_db = os.path.join(output_dir, "kma_db")
    kma_out = os.path.join(output_dir, "kma_out")
    bam_file = os.path.join(output_dir, "mobilome_coverage.bam")
    stats_file = os.path.join(output_dir, "mobilome_coverage.txt")
    
    # Create KMA database and map reads
    run_command(f"kma index -i {regions_fasta} -o {kma_db}", "KMA database creation")
    run_command(
        f"kma -i {reads1} {reads2} -t_db {kma_db} -o {kma_out} "
        f"-sam -mp 20 -mrs 0.75 -bcNano -dense -1t1",
        "KMA read mapping"
    )
    
    # Convert SAM to BAM and sort
    run_command(
        f"samtools view -bS {kma_out}.sam | samtools sort -o {bam_file}",
        "SAM to BAM conversion and sorting"
    )
    run_command(f"samtools index {bam_file}", "BAM indexing")
    
    # Calculate coverage
    run_command(
        f"bedtools genomecov -ibam {bam_file} -d -split > {stats_file}",
        "Coverage calculation"
    )
    
    # Cleanup SAM file
    os.remove(f"{kma_out}.sam")
    
    return bam_file, stats_file

def find_conjunction_reads(bam_file, min_overlap=0.2, max_clips=10, max_distance=1050):
    """Identify conjunction reads from BAM file"""
    logging.info(f"Starting conjunction read analysis on {bam_file}")
    logging.info(f"Parameters: min_overlap={min_overlap}, max_clips={max_clips}, max_distance={max_distance}")
    
    conjunction_reads = set()
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    processed_reads = 0
    qualifying_reads = 0
    
    for read in bam.fetch():
        processed_reads += 1
        if processed_reads % 100000 == 0:
            logging.info(f"Processed {processed_reads} reads, found {qualifying_reads} qualifying reads")
            
        if read.is_paired and not read.is_unmapped:
            overlap = abs(read.template_length) / read.query_length
            clips = len(read.cigartuples) - sum(op[1] for op in read.cigartuples if op[0] == 0)
            
            if (overlap >= min_overlap and 
                clips <= max_clips and 
                abs(read.template_length) <= max_distance):
                conjunction_reads.add(read.query_name)
                qualifying_reads += 1
    
    bam.close()
    logging.info(f"Completed conjunction read analysis:")
    logging.info(f"- Total reads processed: {processed_reads}")
    logging.info(f"- Qualifying reads found: {qualifying_reads}")
    return conjunction_reads

def main():
    parser = argparse.ArgumentParser(description='Complete mobilome analysis')
    parser.add_argument('--contigs', required=True, help='Input contigs FASTA')
    parser.add_argument('--blast-results', required=True, help='BLAST results TSV')
    parser.add_argument('--reads1', required=True, help='Input reads 1')
    parser.add_argument('--reads2', required=True, help='Input reads 2')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--min-coverage', type=float, default=0.7, help='Minimum coverage threshold')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.output_dir)
    
    start_time = datetime.now()
    logger.info(f"Starting mobilome analysis at {start_time}")
    logger.info("Input parameters:")
    for arg, value in vars(args).items():
        logger.info(f"  {arg}: {value}")
    
    try:
        # Extract ARG regions
        logger.info("=== Starting ARG region extraction ===")
        bed_file = os.path.join(args.output_dir, "arg_regions.bed")
        regions_fasta = os.path.join(args.output_dir, "arg_regions.fasta")
        extract_arg_regions(args.contigs, args.blast_results, bed_file, regions_fasta)
        
        # Process mobilome mapping and coverage
        logger.info("=== Starting mobilome mapping and coverage analysis ===")
        bam_file, stats_file = process_mobilome(
            args.reads1, args.reads2, regions_fasta, 
            args.output_dir, args.min_coverage
        )
        
        # Find conjunction reads
        logger.info("=== Starting conjunction read analysis ===")
        conjunction_reads = find_conjunction_reads(bam_file)
        conjunction_file = os.path.join(args.output_dir, "conjunction_reads.txt")
        logger.info(f"Writing {len(conjunction_reads)} conjunction reads to {conjunction_file}")
        with open(conjunction_file, 'w') as f:
            for read in conjunction_reads:
                f.write(f"{read}\n")
        
        # Generate summary
        logger.info("=== Generating final summary ===")
        summary_file = os.path.join(args.output_dir, "mobilome_summary.tsv")
        entries_written = 0
        with open(stats_file) as f_stats, open(summary_file, 'w') as f_out:
            for line in f_stats:
                contig, pos, depth = line.strip().split()
                if float(depth) >= args.min_coverage and contig in conjunction_reads:
                    f_out.write(f"{contig}\t{pos}\t{depth}\n")
                    entries_written += 1
        
        logger.info(f"Wrote {entries_written} entries to summary file")
        
        end_time = datetime.now()
        duration = end_time - start_time
        logger.info(f"Mobilome analysis completed successfully in {duration}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()