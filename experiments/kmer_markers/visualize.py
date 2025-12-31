import sys
import time  # Add missing import
import logging
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
try:
    import seaborn as sns
except ImportError:  # optional
    sns = None
from collections import Counter
try:
    import psutil
except ImportError:  # optional
    psutil = None
from typing import Dict, List, Union, Any
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s.%(msecs)03d [%(levelname)s] %(funcName)s:%(lineno)d - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('visualization.log', mode='a')
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
        logger.debug(f"Starting visualization: {func.__name__}")
        log_memory_usage()
        
        result = func(*args, **kwargs)
        
        end_time = time.time()
        duration = end_time - start_time
        logger.debug(f"Completed visualization: {func.__name__}")
        logger.debug(f"Duration: {duration:.2f} seconds")
        log_memory_usage()
        return result
    return wrapper

@lru_cache(maxsize=128)
def cached_plot_calculation(data_key: str, data: np.ndarray):
    """Cache plot calculations for repeated visualizations."""
    return np.histogram(data, bins='auto')

@log_performance_decorator
def plot_data(data, title="Data Visualization", x_label="X", y_label="Y", output_file=None):
    """
    Basic function to create a plot from input data
    """
    plt.figure(figsize=(10, 6))
    plt.plot(data)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(True)
    if output_file:
        plt.savefig(output_file)
    plt.close()

@log_performance_decorator
def plot_distance_matrix(matrix: np.ndarray, title: str = "Distance Matrix", 
                        output_file: str = None, chunk_size: int = 1000):
    """Optimized distance matrix visualization for large matrices."""
    logger.info(f"Plotting distance matrix with shape: {matrix.shape}")
    
    if matrix.size > chunk_size * chunk_size:
        # For large matrices, downsample or process in chunks
        step = max(1, matrix.shape[0] // chunk_size)
        matrix = matrix[::step, ::step]
    
    # Log matrix statistics
    logger.debug(f"Matrix stats - Min: {matrix.min():.3f}, Max: {matrix.max():.3f}, Mean: {matrix.mean():.3f}")
    
    plt.imshow(matrix)
    plt.colorbar()
    plt.title(title)
    
    if output_file:
        logger.info(f"Saving distance matrix plot to: {output_file}")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.debug(f"Plot saved successfully")
    plt.close()

@log_performance_decorator
def plot_cluster_dendrogram(model, title="Hierarchical Clustering Dendrogram", output_file=None):
    plt.figure(figsize=(10, 7))
    dendrogram(model)
    plt.title(title)
    if output_file:
        plt.savefig(output_file)
    plt.close()

@log_performance_decorator
def plot_motif_logo(motif, title="Motif Logo", output_file=None):
    """Plot a sequence motif logo."""
    plt.figure(figsize=(10, 3))
    
    # Create position-specific probability matrix
    heights = np.array(list(motif.counts.values()))
    heights = heights / heights.sum(axis=0)
    
    # Plot stacked bar chart as simple logo representation
    bottom = np.zeros(len(motif))
    for i, base in enumerate(motif.counts.keys()):
        plt.bar(range(len(motif)), heights[i], bottom=bottom, label=base)
        bottom += heights[i]
    
    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Probability")
    plt.legend()
    
    if output_file:
        plt.savefig(output_file)
    plt.close()

@log_performance_decorator
def plot_marker_statistics(markers, title="Marker Statistics", output_file=None):
    """Plot statistics about the markers."""
    plt.figure(figsize=(12, 6))
    
    # Calculate marker lengths
    lengths = [len(seq) for seq in markers.values()]
    
    # Create histogram
    plt.hist(lengths, bins=30)
    plt.title(title)
    plt.xlabel("Marker Length")
    plt.ylabel("Frequency")
    
    if output_file:
        plt.savefig(output_file)
    plt.close()

@log_performance_decorator
def plot_cluster_sizes(clusters, output_file=None):
    """Plot distribution of cluster sizes."""
    plt.figure(figsize=(10, 6))
    
    # Get cluster sizes
    sizes = [len(members) for members in clusters.values()]
    
    # Create bar plot
    counts = Counter(sizes)
    plt.bar(counts.keys(), counts.values())
    
    plt.title("Cluster Size Distribution")
    plt.xlabel("Cluster Size")
    plt.ylabel("Frequency")
    
    if output_file:
        plt.savefig(output_file)
    plt.close()

@log_performance_decorator
def main():
    """Example usage with performance monitoring."""
    logger.info("Starting visualization example")
    sample_data = np.random.random(100)
    plot_data(sample_data, "Sample Plot", "Time", "Value")
    logger.info("Visualization example completed")

if __name__ == "__main__":
    main()