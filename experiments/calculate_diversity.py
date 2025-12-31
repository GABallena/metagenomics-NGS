import os
import sys
import glob
import pandas as pd
import numpy as np
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.distance import DistanceMatrix
from scipy.stats import entropy
from scipy.optimize import fsolve
import argparse

# Function to read .bracken files and extract species abundance data
def read_bracken(file):
    df = pd.read_csv(file, sep='\t')
    return df['new_est_reads'].values  # Only returning abundance

# Calculate alpha diversity using scikit-bio
def calculate_alpha_diversity(counts, metrics):
    results = {}
    for metric in metrics:
        alpha_values = alpha_diversity(metric, counts)
        results[metric] = alpha_values
    return results

# Custom implementation of Fisher's Alpha
def calculate_fisher_alpha(counts):
    counts = np.array(counts)
    def equation(a):
        return np.sum(np.log((a + np.arange(1, counts.sum()+1)) / a)) - counts.sum()
    
    initial_guess = 1.0
    fisher_alpha = fsolve(equation, initial_guess)[0]
    return fisher_alpha

# Custom implementation of Chao1 (using scipy or numpy)
def calculate_chao1(counts):
    counts = np.array(counts)
    singletons = np.sum(counts == 1)
    doubletons = np.sum(counts == 2)
    if doubletons > 0:
        chao1 = len(counts) + singletons**2 / (2 * doubletons)
    else:
        chao1 = len(counts) + singletons * (singletons - 1) / 2
    return chao1

# Custom implementation of Berger-Parker index
def calculate_berger_parker(counts):
    max_abundance = np.max(counts)
    total_abundance = np.sum(counts)
    return max_abundance / total_abundance

# Calculate beta diversity using scikit-bio
def calculate_beta_diversity(counts, metrics, ids):
    results = {}
    for metric in metrics:
        dm = beta_diversity(metric, counts, ids=ids)
        results[metric] = dm
    return results

def main():
    parser = argparse.ArgumentParser(description="Calculate alpha and beta diversity metrics from .bracken files.")
    parser.add_argument("--input", required=True, help="Directory containing .bracken files.")
    parser.add_argument("--output", required=True, help="Output file for diversity metrics.")
    args = parser.parse_args()

    input_dir = args.input
    output_file = args.output

    # Read the .bracken files
    input_files = glob.glob(os.path.join(input_dir, "*.bracken"))
    if len(input_files) < 2:
        print("Please provide at least two .bracken files for beta diversity calculations.")
        sys.exit(1)

    # Clean sample names
    sample_names = [os.path.basename(f).replace('.bracken','') for f in input_files]
    sample_data = [read_bracken(f) for f in input_files]

    # Pad the samples to the same length
    max_length = max(len(data) for data in sample_data)
    padded_sample_data = [np.pad(data, (0, max_length - len(data)), 'constant') for data in sample_data]
    padded_sample_data = np.array(padded_sample_data)

    # Alpha diversity metrics
    alpha_metrics = ['shannon', 'simpson', 'pielou_e']
    alpha_results = {}
    for i, data in enumerate(padded_sample_data):
        alpha_results[sample_names[i]] = {}
        for metric, value in calculate_alpha_diversity([data], alpha_metrics).items():
            alpha_results[sample_names[i]][metric] = value.values[0]  # Extract the float value
        alpha_results[sample_names[i]]['fisher_alpha'] = calculate_fisher_alpha(data)
        alpha_results[sample_names[i]]['chao1'] = calculate_chao1(data)
        alpha_results[sample_names[i]]['berger_parker'] = calculate_berger_parker(data)

    # Beta diversity metrics
    beta_metrics = ['braycurtis', 'jaccard']
    beta_results = calculate_beta_diversity(padded_sample_data, beta_metrics, sample_names)

    # Write results to the output file
    with open(output_file, "w") as f:
        f.write("Alpha Diversity Results:\n")
        for sample, metrics in alpha_results.items():
            f.write(f"\nSample: {sample}\n")
            for metric, value in metrics.items():
                f.write(f"{metric}: {value}\n")
        f.write("\nBeta Diversity Results:\n")
        for metric, matrix in beta_results.items():
            f.write(f"\n{metric.capitalize()} Distance Matrix:\n")
            distance_matrix = matrix
            f.write(distance_matrix.to_data_frame().to_string() + "\n")

if __name__ == "__main__":
    main()
