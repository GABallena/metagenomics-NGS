import os
import sys
import argparse
import subprocess
from Bio import SeqIO
import math
from concurrent.futures import ProcessPoolExecutor
import logging
import time
try:
    import psutil
except ImportError:  # optional
    psutil = None
from tqdm import tqdm
import multiprocessing as mp
import numpy as np
from functools import lru_cache

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s.%(msecs)03d [%(levelname)s] %(funcName)s:%(lineno)d - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('marker_quantification.log')
    ]
)
logger = logging.getLogger(__name__)

def log_memory_usage():
    """Log current memory usage."""
    if psutil is None:
        return
    process = psutil.Process()
    mem_usage = process.memory_info().rss / 1024 / 1024
    logger.debug(f"Current memory usage: {mem_usage:.2f} MB")

def log_performance_decorator(func):
    """Decorator to log function performance metrics."""
    def wrapper(*args, **kwargs):
        start_time = time.time()
        start_memory = psutil.Process().memory_info().rss
        logger.debug(f"Starting {func.__name__}")
        log_memory_usage()
        
        result = func(*args, **kwargs)
        
        end_time = time.time()
        end_memory = psutil.Process().memory_info().rss
        duration = end_time - start_time
        memory_delta = (end_memory - start_memory) / 1024 / 1024
        
        logger.debug(f"Completed {func.__name__}")
        logger.debug(f"Duration: {duration:.2f} seconds")
        logger.debug(f"Memory delta: {memory_delta:.2f} MB")
        return result
    return wrapper

@log_performance_decorator
def check_dependency(program, flag, name):
    """Ensure a program is installed and accessible."""
    logger.info(f"Checking dependency: {name}")
    try:
        subprocess.run([program, flag], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        logger.info(f"Dependency {name} found and working")
    except subprocess.CalledProcessError:
        logger.error(f"Dependency {name} not found or not working")
        sys.exit(f"Error: {name} is not properly installed or not in PATH.")

@log_performance_decorator
def create_marker_db(marker_file, db_path, program):
    """Create a database for markers using a specified program."""
    logger.info(f"Creating marker database at {db_path} using {program}")
    subprocess.run([program, "-makeudb_usearch", marker_file, "-output", db_path])

@log_performance_decorator
def split_wgs_file(wgs_file, output_dir, reads_per_file):
    """Split large WGS files into smaller files."""
    logger.info(f"Splitting WGS file {wgs_file} into smaller files in {output_dir}")
    file_count = 0
    read_count = 0
    with open(wgs_file, 'r') as infile:
        for line in infile:
            if read_count % reads_per_file == 0:
                if file_count > 0:
                    outfile.close()
                file_count += 1
                outfile = open(f"{output_dir}/split_{file_count}.fasta", "w")
            outfile.write(line)
            read_count += 1
        if not outfile.closed:
            outfile.close()
    logger.info(f"Created {file_count} split files")
    return file_count

@log_performance_decorator
def run_search(marker_db, query_file, output_file, program, threads=1, identity=0.95):
    """Run the search program."""
    logger.info(f"Running search on {query_file} against {marker_db} with {program}")
    cmd = [
        program, "-query", query_file, "-db", marker_db,
        "-out", output_file, "-id", str(identity), "-threads", str(threads)
    ]
    subprocess.run(cmd)

@log_performance_decorator
def calculate_abundance(hit_counts, marker_lengths, total_reads, avg_read_length):
    """Normalize marker hit counts by read length and total reads."""
    logger.info("Calculating abundance of markers")
    abundance = {}
    for marker, count in hit_counts.items():
        marker_length = marker_lengths[marker]
        abundance[marker] = (count / (marker_length * avg_read_length)) * total_reads
    return abundance

@lru_cache(maxsize=1024)
def cached_parse_blast(blast_line: str):
    """Cache blast parsing results."""
    return blast_line.split()[:4]

def process_chunk(chunk_data):
    """Process a chunk of blast output in parallel."""
    hit_counts = {}
    for line in chunk_data:
        query, subject, identity, length = cached_parse_blast(line)
        hit_counts[query] = hit_counts.get(query, 0) + 1
    return hit_counts

@log_performance_decorator
def parse_blast_output(blast_file, chunk_size=10000):
    """Optimized parallel blast output parsing."""
    logger.info(f"Parsing BLAST output from {blast_file}")
    
    # Read file in chunks
    with open(blast_file, 'r') as blast:
        lines = blast.readlines()
    
    # Split into chunks for parallel processing
    chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]
    
    # Process chunks in parallel
    with ProcessPoolExecutor(max_workers=mp.cpu_count()) as executor:
        chunk_results = list(executor.map(process_chunk, chunks))
    
    # Merge results
    hit_counts = {}
    for result in chunk_results:
        for key, value in result.items():
            hit_counts[key] = hit_counts.get(key, 0) + value
    
    return hit_counts

@log_performance_decorator
def process_split_file(args):
    """Process a single split file."""
    db_path, query_file, output_file, program, threads, identity = args
    run_search(db_path, query_file, output_file, program, threads, identity)
    return parse_blast_output(output_file)

@log_performance_decorator
def main():
    # Modified to work with Snakemake
    try:
        start_time = time.time()
        logger.info("Starting marker quantification pipeline")
        log_memory_usage()
        
        if 'snakemake' in globals():
            args = argparse.Namespace(
                markers=snakemake.input.markers,
                output=snakemake.output.quantified,
                program=snakemake.params.program,
                threads=snakemake.params.threads,
                identity=snakemake.config['search_identity']
            )
        else:
            parser = argparse.ArgumentParser(description="Quantify markers using WGS or genome data.")
            parser.add_argument("--markers", required=True, help="Path to marker file.")
            parser.add_argument("--wgs", required=True, help="Path to WGS file.")
            parser.add_argument("--output", required=True, help="Output file for quantified markers.")
            parser.add_argument("--program", default="usearch", help="Program for searching markers.")
            parser.add_argument("--threads", type=int, default=1, help="Number of threads for search.")
            parser.add_argument("--identity", type=float, default=0.95, help="Identity threshold for marker search.")
            args = parser.parse_args()

        # Check dependencies
        check_dependency(args.program, "-h", args.program)

        # Create marker database
        db_path = "markers.udb"
        create_marker_db(args.markers, db_path, args.program)

        # Split WGS file
        output_dir = "tmp_wgs"
        os.makedirs(output_dir, exist_ok=True)
        file_count = split_wgs_file(args.wgs, output_dir, reads_per_file=7000000)

        # Run search on each split file in parallel
        hit_counts = {}
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = []
            for i in range(1, file_count + 1):
                query_file = f"{output_dir}/split_{i}.fasta"
                output_file = f"{output_dir}/results_{i}.out"
                futures.append(executor.submit(process_split_file, (db_path, query_file, output_file, args.program, 1, args.identity)))
            for future in futures:
                file_hits = future.result()
                for key, value in file_hits.items():
                    hit_counts[key] = hit_counts.get(key, 0) + value

        # Normalize abundance
        marker_lengths = {seq.id: len(seq) for seq in SeqIO.parse(args.markers, "fasta")}
        total_reads = sum(hit_counts.values())
        avg_read_length = 100  # Default average read length
        abundance = calculate_abundance(hit_counts, marker_lengths, total_reads, avg_read_length)

        # Write output
        with open(args.output, 'w') as outfile:
            outfile.write("Marker\tAbundance\n")
            for marker, value in abundance.items():
                outfile.write(f"{marker}\t{value}\n")
        
        end_time = time.time()
        total_duration = end_time - start_time
        logger.info(f"Pipeline completed successfully in {total_duration:.2f} seconds")
        log_memory_usage()
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
