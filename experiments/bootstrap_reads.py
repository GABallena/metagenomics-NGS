import os
import random
import subprocess
import sys
import argparse
from glob import glob
import time

# Configuration
DEFAULT_PARAMS = {
    'BOOTSTRAP_FRACTION': 0.1,
    'DEFAULT_SEED': 42,
    'NUM_BOOTSTRAPS': 1,
    'PROGRESS_WIDTH': 50,  # Width of progress bar
}

# Directory structure
DIR_STRUCTURE = {
    'INPUT_PATTERN': "sample_*_*.fq.gz*",
    'BOOTSTRAP_PREFIX': "bootstrap",
    'R1_IDENTIFIER': "_1.fq.gz",
    'R2_IDENTIFIER': "_2.fq.gz"
}
#
# Conda environment
CONDA_ENV = "seqtk_env"
SEQTK_CMD = os.getenv('SEQTK_CMD', 'seqtk')
USE_CONDA = os.getenv('USE_CONDA', '1') == '1'

def run_seqtk_command(cmd, description=""):
    """Run a seqtk command with status reporting."""
    try:
        print(f"Running: {description}...", end=" ", flush=True)
        full_command = f"conda run -n {CONDA_ENV} {cmd}" if USE_CONDA else cmd
        result = subprocess.run(full_command, shell=True, capture_output=True, text=True, check=True)
        print("OK")
        return result
    except subprocess.CalledProcessError as e:
        print("FAIL")
        print(f"Error: {e.stderr}")
        raise

def find_fastq_pairs(input_dir):
    """Find all paired FASTQ files in a directory."""
    r1_files = sorted(glob(os.path.join(input_dir, DIR_STRUCTURE['INPUT_PATTERN'])))
    pairs = []
    seen_samples = set()
    
    for r1 in r1_files:
        # Extract sample name (e.g., "sample_0" from "sample_0_1.fq.gz")
        sample_name = "_".join(os.path.basename(r1).split("_")[:-1])
        if sample_name not in seen_samples:
            r2 = r1.replace(DIR_STRUCTURE['R1_IDENTIFIER'], DIR_STRUCTURE['R2_IDENTIFIER'])
            if os.path.exists(r2):
                pairs.append((r1, r2))
                seen_samples.add(sample_name)
    return pairs

def count_reads(fastq_file):
    """Count total number of reads in a FASTQ file using seqtk."""
    try:
        # Use seqtk fqchk to get stats
        result = run_seqtk_command(
            f"zcat {fastq_file} | wc -l | awk '{{print $1/4}}'",
            f"Counting reads in {os.path.basename(fastq_file)}"
        )
        return int(result.stdout.strip())
    except ValueError as e:
        print(f"Error parsing read count: {e}")
        sys.exit(1)

def bootstrap_reads(r1_file, r2_file, output_r1, output_r2, sample_size, seed=42):
    """Perform bootstrapping by sampling reads with replacement."""
    random.seed(seed)
    
    run_seqtk_command(f"seqtk sample -s{seed} {r1_file} {sample_size} > {output_r1}", f"Bootstrapping {os.path.basename(r1_file)}")
    run_seqtk_command(f"seqtk sample -s{seed} {r2_file} {sample_size} > {output_r2}", f"Bootstrapping {os.path.basename(r2_file)}")

def print_section_header(title, width=None):
    """Print a formatted section header."""
    if width is None:
        width = DEFAULT_PARAMS['PROGRESS_WIDTH']
    print(f"\n{'=' * width}")
    print(f"{title}")
    print(f"{'=' * width}")

def print_progress(message, end="\n"):
    """Print a progress message with timestamp."""
    timestamp = time.strftime('%H:%M:%S')
    print(f"[{timestamp}] {message}", end=end)

def process_bootstrap(pairs, output_dir, bootstrap_fraction, seed, bootstrap_num):
    """Process all pairs for a single bootstrap."""
    print_section_header(f"Processing Bootstrap {bootstrap_num}")
    
    bootstrap_dir = os.path.join(output_dir, f"{DIR_STRUCTURE['BOOTSTRAP_PREFIX']}{bootstrap_num}")
    print_progress(f"Creating output directory: {bootstrap_dir}")
    os.makedirs(bootstrap_dir, exist_ok=True)

    total_processed = 0
    successful = 0
    
    for r1_file, r2_file in pairs:
        print_progress(f"Processing sample: {os.path.basename(r1_file)}")
        
        output_r1 = os.path.join(bootstrap_dir, os.path.basename(r1_file))
        output_r2 = os.path.join(bootstrap_dir, os.path.basename(r2_file))
        
        total_reads = count_reads(r1_file)
        bootstrap_sample_size = int(bootstrap_fraction * total_reads)
        print_progress(f"  Total Reads: {total_reads:,}")
        print_progress(f"  Sample Size: {bootstrap_sample_size:,}")
        
        try:
            bootstrap_reads(r1_file, r2_file, output_r1, output_r2, bootstrap_sample_size, seed)
            total_processed += total_reads
            successful += 1
            
        except Exception as e:
            print(f"Error processing {os.path.basename(r1_file)}: {e}")
            continue
            
    return total_processed, successful

def main():
    parser = argparse.ArgumentParser(
        description="Bootstrap FASTQ reads from a directory",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input-dir", required=True,
                        help="Input directory containing FASTQ files")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Output directory for bootstrapped files")
    parser.add_argument("-f", "--fraction", type=float,
                        default=DEFAULT_PARAMS['BOOTSTRAP_FRACTION'],
                        help="Fraction of reads to sample")
    parser.add_argument("-n", "--num-bootstraps", type=int,
                        default=DEFAULT_PARAMS['NUM_BOOTSTRAPS'],
                        help="Number of bootstrap replicates")
    parser.add_argument("-s", "--seed", type=int,
                        default=DEFAULT_PARAMS['DEFAULT_SEED'],
                        help="Random seed")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    start_time = time.time()
    total_reads_processed = 0
    successful_bootstraps = 0
    failed_bootstraps = 0

    print_section_header("Bootstrap Analysis Configuration")
    print_progress(f"Input directory: {args.input_dir}")
    print_progress(f"Output directory: {args.output_dir}")
    print_progress(f"Bootstrap fraction: {args.fraction}")
    print_progress(f"Number of bootstraps: {args.num_bootstraps}")
    print_progress(f"Random seed: {args.seed}")
    
    print(f"\nStarting bootstrap analysis with {args.num_bootstraps} replicates")
    print("=" * DEFAULT_PARAMS['PROGRESS_WIDTH'])

    # Find paired FASTQ files
    print_section_header("Scanning Input Directory")
    pairs = find_fastq_pairs(args.input_dir)
    if not pairs:
        print_progress("No paired FASTQ files found!")
        sys.exit(1)

    print_progress(f"Found {len(pairs)} paired FASTQ files:")
    for r1, r2 in pairs:
        print_progress(f"  Pair: {os.path.basename(r1)} / {os.path.basename(r2)}")

    # Process each bootstrap
    for bootstrap_num in range(1, args.num_bootstraps + 1):
        start_bootstrap = time.time()
        current_seed = args.seed + bootstrap_num - 1
        
        processed_reads, successful = process_bootstrap(
            pairs, args.output_dir, args.fraction, current_seed, bootstrap_num
        )
        
        total_reads_processed += processed_reads
        successful_bootstraps += successful
        failed_bootstraps += (len(pairs) * 2 - successful)
        
        bootstrap_time = time.time() - start_bootstrap
        print_progress(f"Bootstrap {bootstrap_num} completed in {bootstrap_time:.2f} seconds")

    # Print final statistics
    elapsed_time = time.time() - start_time
    print_section_header("Final Statistics")
    print_progress(f"Total reads processed: {total_reads_processed:,}")
    print_progress(f"Successful bootstraps: {successful_bootstraps}")
    print_progress(f"Failed bootstraps: {failed_bootstraps}")
    print_progress(f"Time elapsed: {elapsed_time:.2f} seconds")
    print_progress(f"Average time per bootstrap: {elapsed_time/args.num_bootstraps:.2f} seconds")

if __name__ == "__main__":
    main()
