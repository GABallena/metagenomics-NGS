import logging
import psutil
import time
import os
import json
import mmap
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Optional, Tuple, List, Any
from scipy.sparse import csr_matrix, vstack, save_npz, load_npz
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
import pystan
import argparse
import signal
import sys
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import normalize
import gc
import psycopg2
from psycopg2.extras import execute_values


# Database configuration (optional)
# For public repos: do not hardcode credentials. Use environment variables instead.
db_config = {
    'dbname': os.getenv('KMER_DB_NAME', 'kmer_db'),
    'user': os.getenv('KMER_DB_USER', 'kmer_user'),
    'password': os.getenv('KMER_DB_PASSWORD', ''),
    'host': os.getenv('KMER_DB_HOST', 'localhost'),
    'port': int(os.getenv('KMER_DB_PORT', '5432')),
}
# Database configuration with connection retries and error handling
def get_db_connection(max_retries=3, retry_delay=2):
    retries = 0
    while retries < max_retries:
        try:
            conn = psycopg2.connect(**db_config)
            conn.set_session(autocommit=True)  # Enable autocommit mode
            return conn
        except psycopg2.OperationalError as e:
            retries += 1
            if retries == max_retries:
                raise Exception(f"Failed to connect to database after {max_retries} attempts: {e}")
            time.sleep(retry_delay)
    return None

# Database connection is optional. If unavailable, fall back to an in-memory vocabulary.
DB_AVAILABLE = False
try:
    conn = get_db_connection()
    if conn is not None:
        with conn.cursor() as cursor:
            cursor.execute("SELECT current_user, current_database()")
            user, db = cursor.fetchone()
            logging.info(f"Connected to database {db} as user {user}")
        conn.close()
        DB_AVAILABLE = True
except Exception as e:
    logging.warning(f"Database connection unavailable; continuing without DB: {e}")
    DB_AVAILABLE = False

# Enhanced logging configuration
logging.basicConfig(
    level=logging.DEBUG,  # Changed from INFO to DEBUG
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('kmer_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

@dataclass
class ProcessingState:
    processed_files: Dict[str, bool]
    current_k: int
    checkpoint_path: str

    def save(self):
        with open(self.checkpoint_path, 'w') as f:
            json.dump({
                'processed_files': self.processed_files,
                'current_k': self.current_k
            }, f)

    @classmethod
    def load(cls, checkpoint_path: str) -> 'ProcessingState':
        if os.path.exists(checkpoint_path):
            with open(checkpoint_path, 'r') as f:
                data = json.load(f)
                return cls(
                    processed_files=data['processed_files'],
                    current_k=data['current_k'],
                    checkpoint_path=checkpoint_path
                )
        return cls({}, 0, checkpoint_path)

def validate_fasta(file_path: Path) -> bool:
    """Validate if a file is a properly formatted FASTA file."""
    try:
        logger.debug(f"Validating FASTA file: {file_path}")
        with open(file_path, 'r') as f:
            for _ in SeqIO.parse(f, 'fasta'):
                pass  # Attempt to parse, raises an exception if invalid
        logger.debug(f"FASTA file is valid: {file_path}")
        return True
    except Exception as e:
        logger.warning(f"Invalid FASTA file skipped: {file_path} ({e})")
        return False

def generate_kmers(sequence: str, k: int) -> List[str]:
    logger.debug(f"Generating k-mers for sequence of length {len(sequence)} with k={k}")
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def create_sparse_row(kmers: List[str], vocab: Dict[str, int]) -> csr_matrix:
    """Create a sparse matrix row for a given list of k-mers and vocabulary."""
    rows, cols, data = [], [], []
    for kmer in kmers:
        if kmer in vocab:
            rows.append(0)
            cols.append(vocab[kmer])
            data.append(1)
    return csr_matrix((data, (rows, cols)), shape=(1, len(vocab)))

def create_sparse_matrix(kmers: List[str], vocab: Dict[str, int]) -> csr_matrix:
    logger.debug(f"Creating sparse matrix for {len(kmers)} k-mers")
    return create_sparse_row(kmers, vocab)

def save_sparse_matrix_to_disk(matrix: csr_matrix, output_path: Path):
    """Save a sparse matrix to disk in .npz format."""
    logger.debug(f"Saving sparse matrix to {output_path}")
    save_npz(output_path, matrix)
    logger.info(f"Sparse matrix saved to {output_path}")

def load_sparse_matrix_from_disk(input_path: Path) -> csr_matrix:
    """Load a sparse matrix from disk."""
    logger.debug(f"Loading sparse matrix from {input_path}")
    matrix = load_npz(input_path)
    logger.info(f"Sparse matrix loaded from {input_path}")
    return matrix

def process_single_file(args: Tuple[str, int, ProcessingState, Path, Dict[str, int]]) -> Optional[Path]:
    try:
        file_path, k_size, state, output_dir, global_vocab = args
        logger.info(f"Processing file: {file_path}")
        if file_path in state.processed_files and state.processed_files[file_path]:
            logger.info(f"Skipping already processed file: {file_path}")
            return None

        output_path = output_dir / f"{Path(file_path).stem}_k{k_size}.npz"
        if output_path.exists():
            logger.info(f"Matrix for {file_path} already exists at {output_path}")
            state.processed_files[file_path] = True
            state.save()
            return output_path

        log_memory_usage()
        with open(file_path, 'r') as f:
            content = f.read()

        kmers = generate_kmers(content, k_size)
        matrix = create_sparse_row(kmers, global_vocab)

        # Save matrix to disk
        save_sparse_matrix_to_disk(matrix, output_path)

        state.processed_files[file_path] = True
        state.save()

        logger.info(f"File processed successfully: {file_path}")
        return output_path

    except Exception as e:
        logger.error(f"Error in worker processing file {file_path}: {e}", exc_info=True)
        return None

def merge_saved_matrices(file_paths: List[Path], global_vocab_size: int, batch_size: int = 100) -> csr_matrix:
    """Merge saved sparse matrices from disk in batches."""
    if not file_paths:
        raise ValueError("No matrices to merge")
    
    if global_vocab_size <= 0:
        raise ValueError(f"Invalid vocabulary size: {global_vocab_size}")
    
    logger.info(f"Merging {len(file_paths)} matrices in batches of {batch_size}")
    result_matrix = None
    
    try:
        for batch_start in range(0, len(file_paths), batch_size):
            batch_end = min(batch_start + batch_size, len(file_paths))
            batch_files = file_paths[batch_start:batch_end]
            
            logger.info(f"Processing batch {batch_start//batch_size + 1}/{(len(file_paths) + batch_size - 1)//batch_size}")
            batch_matrices = []
            
            for fp in batch_files:
                matrix = load_sparse_matrix_from_disk(fp)
                if matrix.shape[1] != global_vocab_size:
                    logger.warning(f"Matrix {fp} has wrong size {matrix.shape[1]}, expected {global_vocab_size}")
                    continue
                batch_matrices.append(matrix)
                
            if batch_matrices:
                batch_combined = vstack(batch_matrices, format='csr')
                if result_matrix is None:
                    result_matrix = batch_combined
                else:
                    result_matrix = vstack([result_matrix, batch_combined], format='csr')
                    
            # Force garbage collection
            batch_matrices = None
            gc.collect()
            logger.debug(f"Memory usage: {psutil.Process().memory_info().rss / 1024 / 1024:.2f} MB")
            
        if result_matrix is None:
            raise ValueError("No valid matrices were processed")
            
        logger.info(f"Final matrix shape: {result_matrix.shape}")
        return result_matrix
        
    except Exception as e:
        logger.error(f"Error merging matrices: {e}")
        raise

def build_global_vocab_inmemory(file_paths: List[str], k_size: int) -> Dict[str, int]:
    """Build a global k-mer vocabulary in memory (public-safe fallback)."""
    vocab = {}
    idx = 0
    for file_path in file_paths:
        try:
            with open(file_path, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    seq = str(record.seq)
                    for kmer in generate_kmers(seq, k_size):
                        if kmer not in vocab:
                            vocab[kmer] = idx
                            idx += 1
        except Exception as e:
            logger.warning(f"Failed to process file for vocab: {file_path} ({e})")
    return vocab

def build_global_vocab_postgres(file_paths: List[str], k_size: int, db_config: Dict[str, str]) -> Dict[str, int]:
    """
    Build a global vocabulary using PostgreSQL with batch inserts.
    Args:
        file_paths: List of input FASTA file paths.
        k_size: K-mer size.
        db_config: PostgreSQL connection configuration.
    Returns:
        Dict[str, int]: A dictionary mapping k-mers to their indices.
    """
    logger.info("Building global vocabulary using PostgreSQL database...")
    conn = psycopg2.connect(**db_config)
    cursor = conn.cursor()

    # Create table if it doesn't exist
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS vocab (
        kmer TEXT PRIMARY KEY,
        kmer_index SERIAL
    );
    """)
    conn.commit()

    # Fix: Properly fetch the count
    cursor.execute("SELECT COUNT(*) FROM vocab")
    current_index = cursor.fetchone()[0]
    logger.info(f"Starting from index {current_index}")

    for file_idx, file_path in enumerate(file_paths, 1):
        logger.info(f"Processing file {file_idx}/{len(file_paths)}: {file_path}")
        try:
            with open(file_path, 'r') as f:
                batch = []
                for record in SeqIO.parse(f, "fasta"):
                    sequence = str(record.seq)
                    kmers = generate_kmers(sequence, k_size)
                    for kmer in kmers:
                        batch.append((kmer,))
                        if len(batch) >= 1000:
                            execute_values(cursor, """
                                INSERT INTO vocab (kmer)
                                VALUES %s
                                ON CONFLICT DO NOTHING
                            """, batch)
                            conn.commit()
                            batch.clear()

                if batch:
                    execute_values(cursor, """
                        INSERT INTO vocab (kmer)
                        VALUES %s
                        ON CONFLICT DO NOTHING
                    """, batch)
                    conn.commit()
                    batch.clear()

        except Exception as e:
            logger.warning(f"Failed to process file {file_path}: {e}")

    logger.info("Finalizing vocabulary...")
    cursor.execute("SELECT kmer, kmer_index FROM vocab")
    vocab_dict = {row[0]: row[1] for row in cursor.fetchall()}

    conn.close()
    logger.info(f"Global vocabulary saved to PostgreSQL database")
    return vocab_dict

def run_refined_bayesian_analysis(matrix):
    """Run refined Bayesian analysis on k-mer data."""
    stan_model_code = """
    data {
        int<lower=0> N; 
        int<lower=0> K; 
        matrix[N, K] X; 
    }
    parameters {
        vector[K] beta;
        real<lower=0> sigma;
    }
    model {
        beta ~ normal(0, 5);
        sigma ~ exponential(1);
        for (n in 1:N)
            X[n, ] ~ normal(beta, sigma);
    }
    """
    stan_data = {
        'N': matrix.shape[0],
        'K': matrix.shape[1],
        'X': normalize(matrix.toarray(), axis=1)
    }
    try:
        sm = pystan.StanModel(model_code=stan_model_code)
        fit = sm.sampling(data=stan_data, iter=2000, warmup=1000, chains=4,
                          control={'adapt_delta': 0.95, 'max_treedepth': 12})
        return fit.extract()['beta'], fit.extract()['sigma']
    except Exception as e:
        logger.error(f"Error in Bayesian analysis: {e}")
        return None, None

def plot_heatmap_of_importance(beta_samples):
    mean_beta = np.mean(beta_samples, axis=0)
    sns.heatmap(mean_beta.reshape(1, -1), cmap="coolwarm", cbar=True)
    plt.title("K-mer Importance Heatmap")
    plt.xlabel("K-mers")
    plt.show()

def plot_clustering(matrix):
    pca = PCA(n_components=2)
    reduced_data = pca.fit_transform(matrix.toarray())
    kmeans = KMeans(n_clusters=3)
    clusters = kmeans.fit_predict(reduced_data)
    plt.scatter(reduced_data[:, 0], reduced_data[:, 1], c=clusters, cmap='viridis')
    plt.title("Sequence Clustering Based on K-mers")
    plt.xlabel("PCA1")
    plt.ylabel("PCA2")
    plt.show()

def plot_kmer_frequency(matrix):
    kmer_counts = matrix.sum(axis=0).A1
    plt.hist(kmer_counts, bins=50, alpha=0.7)
    plt.title("K-mer Frequency Distribution")
    plt.xlabel("K-mer Count")
    plt.ylabel("Frequency")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='K-mer analysis pipeline')
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input directory containing FASTA files')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Output directory')
   
    parser.add_argument('-k', '--kmer_size', type=int, default=21, help='K-mer size (default: 21)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads to use')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    db_path = output_dir / 'vocab.db'
    file_paths = [str(f) for f in input_dir.glob('*.fasta')]

    if not file_paths:
        logger.error("No FASTA files found in input directory")
        return

    logger.info("Building global vocabulary...")
    global_vocab = build_global_vocab_postgres(file_paths, args.kmer_size, db_config) if DB_AVAILABLE else build_global_vocab_inmemory(file_paths, args.kmer_size)

    # Validate vocabulary
    logger.info("Validating vocabulary...")
    if not global_vocab:
        logger.error("Global vocabulary is empty. Check if input FASTA files contain valid sequences.")
        return

    # Validate and debug vocabulary
    vocab_size = len(global_vocab)
    logger.info(f"Vocabulary size: {vocab_size}")
    if vocab_size == 0:
        # Debug: Print first few FASTA files being processed
        for file_path in file_paths[:3]:
            logger.debug(f"Checking file: {file_path}")
            try:
                with open(file_path, 'r') as f:
                    for record in SeqIO.parse(f, "fasta"):
                        logger.debug(f"Found sequence of length {len(record.seq)}")
                        # Generate and print a few sample kmers
                        sample_kmers = generate_kmers(str(record.seq), args.kmer_size)[:5]
                        logger.debug(f"Sample kmers: {sample_kmers}")
                        break
            except Exception as e:
                logger.error(f"Error reading file {file_path}: {e}")
        return

    # Process files in parallel
    logger.info("Processing FASTA files...")
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(process_single_file, (f, args.kmer_size, ProcessingState({}, args.kmer_size, str(checkpoint_dir / 'processing_state.json')), output_dir, global_vocab))
            for f in file_paths
        ]
        for future in as_completed(futures):
            future.result()

    # Merge matrices
    logger.info("Merging sparse matrices...")
    matrix_files = list(output_dir.glob('*.npz'))
    if not matrix_files:
        logger.error("No matrix files found to merge")
        return

    try:
        combined_matrix = merge_saved_matrices(matrix_files, vocab_size, batch_size=100)
    except Exception as e:
        logger.error("Failed to merge matrices", exc_info=True)
        return

    logger.info("Running Bayesian analysis...")
    beta_samples, sigma_samples = run_refined_bayesian_analysis(combined_matrix)

    if beta_samples is not None:
        logger.info("Generating visualizations...")
        plot_heatmap_of_importance(beta_samples)
        plot_clustering(combined_matrix)
        plot_kmer_frequency(combined_matrix)

if __name__ == "__main__":
    main()
