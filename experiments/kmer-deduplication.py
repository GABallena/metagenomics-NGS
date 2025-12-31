import os
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import entropy, f_oneway
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import logging
import yaml
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import shutil
from tqdm import tqdm
import tempfile
import seaborn.objects as so
from matplotlib.gridspec import GridSpec
from statannotations.Annotator import Annotator

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def validate_kmer_sizes(k_values):
    """Validate k-mer sizes."""
    for k in k_values:
        if not isinstance(k, int) or k < 1:
            raise ValueError(f"Invalid k-mer size: {k}. Must be positive integer.")
        if k > 100:  # practical limit for most applications
            logging.warning(f"Large k-mer size ({k}) may require significant memory")

def cleanup_temp_files(temp_files):
    """Clean up temporary files."""
    for file in temp_files:
        try:
            if os.path.exists(file):
                os.remove(file)
        except Exception as e:
            logging.warning(f"Failed to remove temporary file {file}: {e}")

# Enhanced Jellyfish runner with memory management
def run_jellyfish(input_file, output_file, k, temp_dir):
    try:
        with tempfile.NamedTemporaryFile(dir=temp_dir, suffix='.fasta', delete=False) as temp_fasta:
            fasta_file = temp_fasta.name
            subprocess.run(["seqtk", "seq", "-A", input_file, "-o", fasta_file], 
                         check=True, capture_output=True, text=True)

        # Estimate memory requirements
        file_size = os.path.getsize(input_file)
        mem_estimate = max(int(file_size * 1.5 / 1e6), 100)  # minimum 100M
        
        subprocess.run([
            "jellyfish", "count",
            "-m", str(k),
            "-s", f"{mem_estimate}M",
            "-t", str(os.cpu_count() or 1),
            fasta_file, "-o", output_file
        ], check=True, capture_output=True, text=True)

        dump_file = output_file + ".dump"
        subprocess.run(["jellyfish", "dump", "-c", output_file, "-o", dump_file],
                     check=True, capture_output=True, text=True)
        
        return dump_file, [fasta_file, output_file]
    except subprocess.CalledProcessError as e:
        logging.error(f"Jellyfish error: {e.stderr}")
        raise
    except Exception as e:
        logging.error(f"Error in Jellyfish pipeline: {e}")
        raise

# Function to calculate Shannon entropy of k-mers
def calculate_entropy(dump_file):
    """Calculate Shannon entropy of k-mers with proper type handling."""
    try:
        kmer_counts = pd.read_csv(dump_file, sep=" ", header=None, names=["kmer", "count"])
        if kmer_counts.empty:
            logging.warning(f"Empty k-mer counts file: {dump_file}")
            return 0.0
            
        # Convert counts to float64 and ensure positive values
        counts = kmer_counts["count"].astype(np.float64)
        total_count = counts.sum()
        
        if total_count == 0:
            logging.warning(f"Zero total counts in {dump_file}")
            return 0.0
            
        # Calculate probabilities and ensure they sum to 1
        probabilities = counts / total_count
        probabilities = probabilities[probabilities > 0]  # Remove zero probabilities
        
        if len(probabilities) == 0:
            logging.warning(f"No valid probabilities in {dump_file}")
            return 0.0
            
        # Calculate entropy manually if scipy fails
        try:
            return float(entropy(probabilities))
        except:
            # Manual entropy calculation as fallback
            return float(-np.sum(probabilities * np.log2(probabilities)))
            
    except Exception as e:
        logging.error(f"Error calculating entropy: {e}")
        raise

# Function to calculate additional metrics
def calculate_additional_metrics(dump_file):
    try:
        kmer_counts = pd.read_csv(dump_file, sep=" ", header=None, names=["kmer", "count"])
        total_count = kmer_counts["count"].sum()
        unique_kmers = len(kmer_counts)
        max_kmer_count = kmer_counts["count"].max()
        mean_kmer_count = kmer_counts["count"].mean()
        
        # Calculate Simpson's Index
        probabilities = kmer_counts["count"] / total_count
        simpson_index = 1 - sum(probabilities**2)
        
        # Calculate Chao1 estimator
        singletons = sum(kmer_counts["count"] == 1)
        doubletons = sum(kmer_counts["count"] == 2)
        if doubletons > 0:
            chao1 = unique_kmers + (singletons**2 / (2 * doubletons))
        else:
            chao1 = unique_kmers + singletons
        
        return {
            "Unique_kmers": unique_kmers,
            "Max_kmer_count": max_kmer_count,
            "Mean_kmer_count": mean_kmer_count,
            "Total_kmers": total_count,
            "Simpson_index": simpson_index,
            "Chao1_estimator": chao1
        }
    except Exception as e:
        logging.error(f"Error calculating additional metrics: {e}")
        raise

def process_single_file(args):
    """Process a single file for parallel execution."""
    input_file, k, dataset_name, temp_dir = args
    try:
        output_file = os.path.join(temp_dir, f"{dataset_name}_k{k}.jf")
        dump_file, temp_files = run_jellyfish(input_file, output_file, k, temp_dir)
        H = calculate_entropy(dump_file)
        metrics = calculate_additional_metrics(dump_file)
        cleanup_temp_files(temp_files + [dump_file])
        return {
            "Dataset": dataset_name,
            "k": k,
            "Entropy": H,
            **metrics
        }
    except Exception as e:
        logging.error(f"Failed to process {dataset_name} with k={k}: {e}")
        return None

# Enhanced entropy analysis with parallel processing
def entropy_analysis(input_dirs, k_values):
    results = []
    temp_dir = tempfile.mkdtemp()
    
    try:
        # Prepare arguments for parallel processing
        process_args = []
        for dataset_type, folder in input_dirs.items():
            if not os.path.isdir(folder):
                logging.error(f"Directory not found: {folder}")
                continue
            for file_name in os.listdir(folder):
                if file_name.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                    input_file = os.path.join(folder, file_name)
                    dataset_name = f"{dataset_type}_{Path(file_name).stem}"
                    for k in k_values:
                        process_args.append((input_file, k, dataset_name, temp_dir))

        # Process files in parallel
        with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            for result in tqdm(executor.map(process_single_file, process_args),
                             total=len(process_args),
                             desc="Processing files"):
                if result:
                    results.append(result)

    finally:
        # Cleanup temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    return pd.DataFrame(results)

def setup_plotting_style():
    """Configure global plotting style"""
    plt.style.use('seaborn-v0_8-paper')
    sns.set_palette("husl")
    plt.rcParams.update({
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'font.size': 10,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10
    })

def plot_entropy(df, output_dir):
    """Enhanced entropy distribution visualization"""
    plt.figure(figsize=(12, 8))
    
    # Create main entropy plot
    ax = sns.boxplot(data=df, x="k", y="Entropy", hue="Dataset", 
                    palette="husl", showfliers=False)
    
    # Add individual points with jitter
    sns.stripplot(data=df, x="k", y="Entropy", hue="Dataset", 
                 dodge=True, alpha=0.3, size=4)
    
    # Add statistical annotations
    pairs = [((k, "deduplicated"), (k, "non_deduplicated")) for k in df['k'].unique()]
    annotator = Annotator(ax, pairs, data=df, x="k", y="Entropy", hue="Dataset")
    annotator.configure(test='t-test_ind', text_format='star', loc='outside')
    annotator.apply_and_annotate()
    
    plt.title("Shannon Entropy Distribution Across k-mer Sizes", 
              fontsize=16, pad=20)
    plt.xlabel("k-mer Size", fontsize=14)
    plt.ylabel("Shannon Entropy", fontsize=14)
    plt.legend(title="Dataset", title_fontsize=12, fontsize=10)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/entropy_distribution.png", bbox_inches='tight')
    plt.close()

def plot_metrics(df, output_dir):
    """Enhanced metrics visualization with subplots"""
    metrics = {
        "Diversity Metrics": ["Simpson_index", "Chao1_estimator"],
        "Count Metrics": ["Unique_kmers", "Total_kmers"],
        "Statistical Metrics": ["Mean_kmer_count", "Max_kmer_count"]
    }
    
    for group_name, metric_group in metrics.items():
        fig = plt.figure(figsize=(15, 8))
        gs = GridSpec(2, 2, figure=fig)
        
        # Main plots
        for i, metric in enumerate(metric_group):
            ax = fig.add_subplot(gs[i])
            
            # Create violin plot with boxplot inside
            sns.violinplot(data=df, x="k", y=metric, hue="Dataset",
                         split=True, inner="box", ax=ax)
            
            # Add individual points
            sns.stripplot(data=df, x="k", y=metric, hue="Dataset",
                        dodge=True, alpha=0.3, size=3, ax=ax)
            
            ax.set_title(f"{metric.replace('_', ' ').title()}", fontsize=12)
            ax.set_xlabel("k-mer Size", fontsize=10)
            ax.set_ylabel(metric.replace('_', ' ').title(), fontsize=10)
            
            # Only show legend for first subplot
            if i != 0:
                ax.get_legend().remove()
        
        # Add correlation heatmap
        ax_corr = fig.add_subplot(gs[:, -1])
        correlation_data = df[metric_group + ['k']].corr()
        sns.heatmap(correlation_data, annot=True, cmap='coolwarm', 
                   center=0, ax=ax_corr)
        ax_corr.set_title("Correlation Matrix", fontsize=12)
        
        plt.suptitle(f"{group_name} Analysis", fontsize=16, y=1.02)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{group_name.lower().replace(' ', '_')}.png",
                   bbox_inches='tight', dpi=300)
        plt.close()

def plot_summary_dashboard(df, output_dir):
    """Create a summary dashboard of all metrics"""
    fig = plt.figure(figsize=(20, 12))
    gs = GridSpec(3, 3, figure=fig)
    
    metrics = ["Entropy", "Simpson_index", "Chao1_estimator", 
               "Unique_kmers", "Total_kmers", "Mean_kmer_count"]
    
    for i, metric in enumerate(metrics):
        ax = fig.add_subplot(gs[i // 3, i % 3])
        
        sns.boxplot(data=df, x="k", y=metric, hue="Dataset", 
                   showfliers=False, ax=ax)
        
        ax.set_title(metric.replace('_', ' ').title(), fontsize=12)
        ax.set_xlabel("k-mer Size", fontsize=10)
        ax.tick_params(labelsize=8)
        
        if i % 3 != 0:
            ax.get_legend().remove()
    
    plt.suptitle("K-mer Analysis Dashboard", fontsize=16, y=0.95)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/summary_dashboard.png", 
                bbox_inches='tight', dpi=300)
    plt.close()

def save_to_yaml(results_df, yaml_path):
    """Save analysis results to YAML format."""
    results_dict = {
        'summary': {
            'total_samples': len(results_df),
            'k_values': sorted(results_df['k'].unique().tolist()),
            'datasets': sorted(results_df['Dataset'].unique().tolist())
        },
        'metrics': results_df.to_dict(orient='records')
    }
    
    with open(yaml_path, 'w') as f:
        yaml.dump(results_dict, f, default_flow_style=False)
    logging.info(f"Results saved to YAML: {yaml_path}")

# Main script
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze k-mer entropy and additional metrics.")
    parser.add_argument("-d", required=True, help="Path to deduplicated FASTQ files.")
    parser.add_argument("-nd", required=True, help="Path to non-deduplicated FASTQ files.")
    parser.add_argument("-k", type=int, nargs="+", default=[21, 31, 41],
                       help="List of k-mer sizes to analyze")
    parser.add_argument("-c", "--config", help="Path to config file")
    parser.add_argument("--output-dir", default="results",
                       help="Directory for output files")
    parser.add_argument("-y", "--yaml", help="Path to save results in YAML format")
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Setup logging
    log_file = os.path.join(args.output_dir, "analysis.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

    # Load config if provided
    if args.config:
        with open(args.config) as f:
            config = yaml.safe_load(f)
            k_values = config.get('k_values', args.k)
    else:
        k_values = args.k

    # Validate k-mer sizes
    validate_kmer_sizes(k_values)

    # Directories containing datasets
    input_dirs = {
        "deduplicated": args.d,
        "non_deduplicated": args.nd
    }

    # Validate input directories
    for key, path in input_dirs.items():
        if not os.path.isdir(path):
            logging.error(f"Invalid directory for {key}: {path}")
            exit(1)

    # Ensure required tools are installed
    for tool in ["jellyfish", "seqtk"]:
        if subprocess.run(["which", tool], capture_output=True).returncode != 0:
            logging.error(f"{tool} is not installed or not in PATH.")
            exit(1)

    # Perform entropy analysis
    entropy_df = entropy_analysis(input_dirs, k_values)

    # Save results to a CSV file
    entropy_df.to_csv("kmer_entropy_results.csv", index=False)

    # Save results to YAML if specified
    if args.yaml:
        save_to_yaml(entropy_df, args.yaml)

    # Perform ANOVA
    for k in k_values:
        subset = entropy_df[entropy_df["k"] == k]
        f_stat, p_val = f_oneway(
            subset[subset["Dataset"].str.contains("deduplicated")]["Entropy"],
            subset[subset["Dataset"].str.contains("non_deduplicated")]["Entropy"]
        )
        logging.info(f"ANOVA for k={k}: F={f_stat:.2f}, p={p_val:.3e}")

    # Setup plotting style
    setup_plotting_style()
    
    # Create visualizations
    plot_entropy(entropy_df, args.output_dir)
    plot_metrics(entropy_df, args.output_dir)
    plot_summary_dashboard(entropy_df, args.output_dir)
