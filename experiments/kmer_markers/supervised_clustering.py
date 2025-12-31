import numpy as np
from Bio import SeqIO
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances, silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from tqdm import tqdm
import logging
import time
from datetime import datetime
import os
import sys
import json
import hashlib
import warnings
from typing import List, Tuple, Dict, Any
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
import pickle
from pathlib import Path
from logging.handlers import RotatingFileHandler
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import squareform

@dataclass
class ClusteringConfig:
    """Configuration parameters with validation."""
    fasta_file: str
    similarity_threshold: float
    num_clusters: int
    dbscan_eps: float
    dbscan_min_samples: int
    test_size: float

    @classmethod
    def from_dict(cls, config_dict: dict) -> 'ClusteringConfig':
        """Create instance from dictionary, ignoring unknown parameters."""
        valid_params = cls.__dataclass_fields__.keys()
        filtered_dict = {k: v for k, v in config_dict.items() if k in valid_params}
        return cls(**filtered_dict)
    
    def validate(self) -> None:
        """Validate configuration parameters."""
        if not os.path.exists(self.fasta_file):
            raise FileNotFoundError(f"FASTA file not found: {self.fasta_file}")
        if not (0 < self.similarity_threshold <= 1):
            raise ValueError("Similarity threshold must be between 0 and 1")
        if self.num_clusters < 2:
            raise ValueError("Number of clusters must be at least 2")
        if self.dbscan_eps <= 0:
            raise ValueError("DBSCAN epsilon must be positive")
        if self.dbscan_min_samples < 1:
            raise ValueError("DBSCAN min_samples must be at least 1")
        if not (0 < self.test_size < 1):
            raise ValueError("Test size must be between 0 and 1")

def compute_similarity(seq1: str, seq2: str) -> float:
    """
    Compute the similarity between two sequences.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        
    Returns:
        float: Similarity score between 0 and 1
    """
    if not seq1 or not seq2:
        return 0.0
    
    try:
        matches = sum(a == b for a, b in zip(seq1, seq2))
        total = min(len(seq1), len(seq2))
        return matches / total if total > 0 else 0.0
    except Exception as e:
        logging.error(f"Error computing similarity: {str(e)}")
        return 0.0

class RobustClusteringPipeline:
    def __init__(self, config: ClusteringConfig):
        self.config = config
        self.setup_logging()
        self.setup_backup_dir()
        self.sequences = []
        self.sequence_lengths = []  # Add storage for sequences
        
    def setup_logging(self) -> None:
        """Configure robust logging with rotation and backup."""
        try:
            log_dir = Path("logs")
            log_dir.mkdir(exist_ok=True)
            
            # Configure root logger
            root_logger = logging.getLogger()
            root_logger.setLevel(logging.INFO)
            
            # Create formatters and handlers
            formatter = logging.Formatter(
                '%(asctime)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'
            )
            
            # File handler with rotation
            file_handler = RotatingFileHandler(
                log_dir / "clustering_analysis.log",
                maxBytes=10_000_000,  # 10MB
                backupCount=5
            )
            file_handler.setFormatter(formatter)
            
            # Console handler
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            
            # Remove existing handlers if any
            root_logger.handlers.clear()
            
            # Add handlers
            root_logger.addHandler(file_handler)
            root_logger.addHandler(console_handler)
            
            # Log uncaught exceptions
            def exception_handler(exc_type, exc_value, exc_traceback):
                root_logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
            sys.excepthook = exception_handler
            
            root_logger.info("Logging system initialized")
            
        except Exception as e:
            print(f"Error setting up logging: {str(e)}")
            raise

    def setup_backup_dir(self) -> None:
        """Create backup directory for intermediate results."""
        self.backup_dir = Path("backup")
        self.backup_dir.mkdir(exist_ok=True)

    def compute_checksum(self, data: Any) -> str:
        """Compute checksum for data validation."""
        return hashlib.md5(str(data).encode()).hexdigest()

    def save_checkpoint(self, data: Any, name: str) -> None:
        """Save intermediate results with checksum."""
        checkpoint_path = self.backup_dir / f"{name}.pkl"
        checksum = self.compute_checksum(data)
        
        with open(checkpoint_path, 'wb') as f:
            pickle.dump({'data': data, 'checksum': checksum}, f)
        
        logging.info(f"Checkpoint saved: {name} (checksum: {checksum})")

    def load_checkpoint(self, name: str) -> Any:
        """Load checkpoint with validation."""
        checkpoint_path = self.backup_dir / f"{name}.pkl"
        
        if checkpoint_path.exists():
            with open(checkpoint_path, 'rb') as f:
                checkpoint = pickle.load(f)
                
            stored_checksum = checkpoint['checksum']
            computed_checksum = self.compute_checksum(checkpoint['data'])
            
            if stored_checksum != computed_checksum:
                raise ValueError(f"Checkpoint validation failed for {name}")
                
            return checkpoint['data']
        return None

    def process_sequences(self) -> Tuple[List[str], List[str], np.ndarray]:
        """Process sequences with validation and error handling."""
        try:
            sequences = []
            sequence_ids = []
            sequence_lengths = []  # Add local variable for lengths
            
            with open(self.config.fasta_file, 'r') as f:
                content = f.read()
                if not content.strip():
                    raise ValueError("Empty FASTA file")
                    
            def normalize_sequence(seq: str) -> str:
                """Normalize sequence by converting to uppercase and handling special characters."""
                # Convert to uppercase
                seq = seq.upper()
                # Replace common ambiguous nucleotides with N
                ambiguous_map = {
                    'R': 'N',  # A or G
                    'Y': 'N',  # C or T
                    'S': 'N',  # G or C
                    'W': 'N',  # A or T
                    'K': 'N',  # G or T
                    'M': 'N',  # A or C
                    'B': 'N',  # C or G or T
                    'D': 'N',  # A or G or T
                    'H': 'N',  # A or C or T
                    'V': 'N',  # A or C or G
                }
                for k, v in ambiguous_map.items():
                    seq = seq.replace(k, v)
                # Replace any remaining non-standard characters with N
                return ''.join(c if c in {'A', 'T', 'G', 'C', 'N', '-'} else 'N' for c in seq)

            for record in SeqIO.parse(self.config.fasta_file, "fasta"):
                seq = str(record.seq)
                if not seq:
                    logging.warning(f"Empty sequence found for {record.id}")
                    continue
                
                # Normalize sequence
                normalized_seq = normalize_sequence(seq)
                if normalized_seq.count('N') / len(normalized_seq) > 0.5:
                    logging.warning(f"Sequence {record.id} has >50% ambiguous bases")
                    continue
                    
                sequences.append(normalized_seq)
                sequence_ids.append(record.id)

            if not sequences:
                raise ValueError("No valid sequences found")
                
            # Store sequences and lengths in instance variables
            self.sequences = sequences
            self.sequence_lengths = [len(s) for s in sequences]
                
            logging.info(f"Processed {len(sequences)} sequences after filtering")
            return sequences, sequence_ids, np.array(self.sequence_lengths)
            
        except Exception as e:
            logging.error(f"Error processing sequences: {str(e)}")
            raise

    def compute_similarity_matrix(self, sequences: List[str]) -> np.ndarray:
        """Compute similarity matrix with parallel processing and validation."""
        try:
            n = len(sequences)
            similarity_matrix = np.zeros((n, n))
            
            def compute_batch(batch: List[Tuple[int, int]]) -> List[Tuple[int, int, float]]:
                return [(i, j, compute_similarity(sequences[i], sequences[j])) 
                        for i, j in batch]

            # Generate batch indices
            batches = [
                [(i, j) for j in range(i, n)]
                for i in range(n)
            ]
            
            # Parallel processing
            with ThreadPoolExecutor() as executor:
                results = list(tqdm(
                    executor.map(compute_batch, batches),
                    total=n,
                    desc="Computing similarities"
                ))

            # Fill matrix
            for batch_result in results:
                for i, j, sim in batch_result:
                    similarity_matrix[i, j] = sim
                    similarity_matrix[j, i] = sim

            # Validate matrix
            if not np.allclose(similarity_matrix, similarity_matrix.T):
                raise ValueError("Similarity matrix is not symmetric")
            if not np.allclose(np.diag(similarity_matrix), 1.0):
                raise ValueError("Diagonal elements are not 1.0")

            return similarity_matrix

        except Exception as e:
            logging.error(f"Error computing similarity matrix: {str(e)}")
            raise

    def run_clustering(self) -> Dict[str, Any]:
        """Run clustering pipeline with comprehensive error handling."""
        try:
            self.config.validate()
            results = {}
            
            # Process sequences
            sequences, sequence_ids, lengths = self.process_sequences()
            self.save_checkpoint((sequences, sequence_ids), "sequences")
            
            # Compute similarity matrix
            similarity_matrix = self.compute_similarity_matrix(sequences)
            self.save_checkpoint(similarity_matrix, "similarity_matrix")
            
            # Run clustering algorithms with validation
            clustering_methods = {
                'cast': lambda: self.run_cast_clustering(similarity_matrix),
                'kmeans': lambda: self.run_kmeans_clustering(similarity_matrix),
                'hierarchical': lambda: self.run_hierarchical_clustering(similarity_matrix),
                'dbscan': lambda: self.run_dbscan_clustering(similarity_matrix)
            }
            
            for method, func in clustering_methods.items():
                try:
                    results[method] = func()
                    self.save_checkpoint(results[method], f"results_{method}")
                except Exception as e:
                    logging.error(f"Error in {method} clustering: {str(e)}")
                    results[method] = None

            # Validate results
            self.validate_results(results)
            
            return results
            
        except Exception as e:
            logging.error(f"Pipeline failed: {str(e)}")
            raise

    def run_dbscan_clustering(self, similarity_matrix: np.ndarray) -> Dict[str, Any]:
        """Perform DBSCAN clustering."""
        try:
            logging.info("Starting DBSCAN clustering")
            start_time = time.time()
            
            dbscan = DBSCAN(eps=self.config.dbscan_eps, min_samples=self.config.dbscan_min_samples)
            labels = dbscan.fit_predict(1 - similarity_matrix)  # Use distance matrix
            
            results = {
                'labels': labels.tolist(),
                'n_clusters': len(set(labels)) - (1 if -1 in labels else 0),
                'n_noise': list(labels).count(-1)
            }
            
            if len(set(labels)) > 1:
                results['silhouette'] = float(silhouette_score(similarity_matrix, labels))
            
            logging.info(f"DBSCAN completed in {time.time() - start_time:.2f} seconds")
            return results
            
        except Exception as e:
            logging.error(f"Error in DBSCAN clustering: {str(e)}")
            raise

    def run_cast_clustering(self, similarity_matrix: np.ndarray) -> Dict[str, Any]:
        """Perform CAST clustering using K-means approximation."""
        try:
            logging.info("Starting CAST clustering")
            start_time = time.time()
            
            kmeans = KMeans(
                n_clusters=self.config.num_clusters, 
                random_state=42,
                n_init=10
            )
            labels = kmeans.fit_predict(similarity_matrix)
            
            results = {
                'labels': labels.tolist(),
                'n_clusters': len(set(labels)),
                'cluster_centers': kmeans.cluster_centers_.tolist(),
                'inertia': float(kmeans.inertia_)
            }
            
            if len(set(labels)) > 1:
                results['silhouette'] = float(silhouette_score(similarity_matrix, labels))
            
            logging.info(f"CAST clustering completed in {time.time() - start_time:.2f} seconds")
            return results
            
        except Exception as e:
            logging.error(f"Error in CAST clustering: {str(e)}")
            raise

    def run_kmeans_clustering(self, similarity_matrix: np.ndarray) -> Dict[str, Any]:
        """Perform K-means clustering."""
        try:
            logging.info("Starting K-means clustering")
            start_time = time.time()
            
            kmeans = KMeans(
                n_clusters=self.config.num_clusters,
                random_state=42,
                n_init=10,
                max_iter=300
            )
            labels = kmeans.fit_predict(similarity_matrix)
            
            results = {
                'labels': labels.tolist(),
                'n_clusters': len(set(labels)),
                'cluster_centers': kmeans.cluster_centers_.tolist(),
                'inertia': float(kmeans.inertia_)
            }
            
            if len(set(labels)) > 1:
                results['silhouette'] = float(silhouette_score(similarity_matrix, labels))
            
            logging.info(f"K-means clustering completed in {time.time() - start_time:.2f} seconds")
            return results
            
        except Exception as e:
            logging.error(f"Error in K-means clustering: {str(e)}")
            raise

    def run_hierarchical_clustering(self, similarity_matrix: np.ndarray) -> Dict[str, Any]:
        """Perform hierarchical clustering."""
        try:
            logging.info("Starting hierarchical clustering")
            start_time = time.time()
            
            # Convert similarity to distance and condense matrix
            distance_matrix = 1 - similarity_matrix
            condensed_dist = squareform(distance_matrix)  # Add this import from scipy.spatial.distance
            
            # Compute linkage matrix
            linkage_matrix = linkage(condensed_dist, method='average')
            
            # Cut tree to get clusters
            labels = fcluster(
                linkage_matrix,
                t=self.config.num_clusters,
                criterion='maxclust'
            )
            
            results = {
                'labels': labels.tolist(),
                'n_clusters': len(set(labels)),
                'linkage_matrix': linkage_matrix.tolist()
            }
            
            if len(set(labels)) > 1:
                results['silhouette'] = float(silhouette_score(similarity_matrix, labels))
            
            logging.info(f"Hierarchical clustering completed in {time.time() - start_time:.2f} seconds")
            return results
            
        except Exception as e:
            logging.error(f"Error in hierarchical clustering: {str(e)}")
            raise

    def validate_results(self, results: Dict[str, Any]) -> None:
        """Validate clustering results."""
        try:
            if not results:
                raise ValueError("No clustering results available")
                
            # Check if at least one clustering method succeeded
            successful_methods = [method for method, result in results.items() if result is not None]
            if not successful_methods:
                raise ValueError("All clustering methods failed")
                
            logging.info(f"Successful clustering methods: {successful_methods}")
            
        except Exception as e:
            logging.error(f"Error validating results: {str(e)}")
            raise

DEFAULT_CONFIG = {
    "fasta_file": "data/input_sequences.fasta",
    "similarity_threshold": 0.9,
    "num_clusters": 3,
    "dbscan_eps": 0.1,
    "dbscan_min_samples": 2,
    "test_size": 0.2
}

if __name__ == "__main__":
    try:
        # Try to load configuration from file, fall back to defaults if not found
        config_dict = DEFAULT_CONFIG.copy()  # Make a copy of defaults
        config_path = Path("config.json")
        
        if config_path.exists():
            try:
                with open(config_path) as f:
                    file_config = json.load(f)
                    # Only update with relevant configuration keys
                    for key in DEFAULT_CONFIG.keys():
                        if key in file_config:
                            config_dict[key] = file_config[key]
                    logging.info("Loaded configuration from config.json")
            except json.JSONDecodeError as e:
                logging.warning(f"Error reading config.json: {e}. Using default configuration.")
        else:
            logging.warning("config.json not found. Using default configuration.")
            # Create default config file
            with open(config_path, 'w') as f:
                json.dump(DEFAULT_CONFIG, f, indent=4)
            logging.info("Created default config.json file")
        
        # Create configuration object
        try:
            config = ClusteringConfig(**config_dict)
            config.validate()  # Validate configuration
        except Exception as e:
            logging.error(f"Invalid configuration: {str(e)}")
            raise
            
        pipeline = RobustClusteringPipeline(config)
        results = pipeline.run_clustering()
        
        # Save final results
        results_output = {
            'clustering_results': results,
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'config': config.__dict__,
                'n_sequences': len(pipeline.sequences),
                'avg_sequence_length': float(np.mean(pipeline.sequence_lengths))
            }
        }
        
        with open(os.path.join(config.output_dir, "results.json"), "w") as f:
            json.dump(results_output, f, indent=2)
        
        logging.info(f"Results saved to {os.path.join(config.output_dir, 'results.json')}")
            
    except Exception as e:
        logging.critical(f"Fatal error: {str(e)}")
        sys.exit(1)

