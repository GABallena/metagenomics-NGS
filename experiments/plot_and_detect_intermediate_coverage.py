import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Function to calculate GC content
def calculate_gc(sequence):
    g = sequence.count('G')
    c = sequence.count('C')
    return (g + c) / len(sequence)

# Function to calculate tetranucleotide frequencies
def calculate_tetra_freq(sequence):
    tetranucleotides = ['AAAA', 'AAAC', 'AAAG', '...']  # Add all tetranucleotides
    freq = {tetra: 0 for tetra in tetranucleotides}
    
    for i in range(len(sequence) - 3):
        tetra = sequence[i:i+4]
        if tetra in freq:
            freq[tetra] += 1
    return np.array(list(freq.values())) / len(sequence)

# Function to calculate codon usage bias
def calculate_codon_bias(sequence):
    codons = ['ATG', 'AAA', 'AAC', '...']  # Add all codons
    bias = {codon: 0 for codon in codons}
    
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in bias:
            bias[codon] += 1
    return np.array(list(bias.values())) / len(sequence)

# Sliding window function
def sliding_window(sequence, window_size, step_size):
    for i in range(0, len(sequence) - window_size + 1, step_size):
        yield sequence[i:i + window_size]

# Main function to analyze each contig
def analyze_contig(contig, window_size, step_size, output_plot):
    gc_content = []
    tetra_freqs = []
    codon_bias = []
    
    for window in sliding_window(contig, window_size, step_size):
        # Calculate GC content, tetranucleotide frequencies, and codon usage bias
        gc_content.append(calculate_gc(window))
        tetra_freqs.append(calculate_tetra_freq(window))
        codon_bias.append(calculate_codon_bias(window))
    
    # Convert lists to numpy arrays
    gc_content = np.array(gc_content)
    tetra_freqs = np.array(tetra_freqs)
    codon_bias = np.array(codon_bias)
    
    # Perform one-way ANOVA on GC content, tetranucleotide frequencies, and codon usage bias
    anova_gc = stats.f_oneway(*[gc_content[i:i + 3] for i in range(0, len(gc_content), 3)])
    anova_tetra = stats.f_oneway(*[tetra_freqs[i:i + 3].mean(axis=0) for i in range(0, len(tetra_freqs), 3)])
    anova_codon = stats.f_oneway(*[codon_bias[i:i + 3].mean(axis=0) for i in range(0, len(codon_bias), 3)])
    
    print(f"ANOVA GC: {anova_gc.pvalue}, ANOVA Tetra: {anova_tetra.pvalue}, ANOVA Codon: {anova_codon.pvalue}")
    
    # Combine all features for PCA
    combined_features = np.column_stack([gc_content, tetra_freqs.mean(axis=1), codon_bias.mean(axis=1)])
    
    # PCA analysis
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(combined_features)
    
    # Plot PCA results
    plt.figure(figsize=(10, 6))
    plt.scatter(principal_components[:, 0], principal_components[:, 1], c='blue', label='PCA')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA of Sliding Window Features')
    plt.savefig(output_plot)

if __name__ == "__main__":
    # Input contig and output file
    contig_file = sys.argv[1]
    output_plot = sys.argv[2]
    
    # Parameters
    window_size = 1000  # You can adjust
    step_size = 500     # You can adjust
    
    # Read contig sequence from file
    with open(contig_file) as f:
        contig = f.read().strip()
    
    # Analyze the contig with sliding window, ANOVA, and PCA
    analyze_contig(contig, window_size, step_size, output_plot)
