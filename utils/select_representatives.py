#!/usr/bin/env python3
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from typing import Dict, Set, List, Optional
import subprocess
from Bio.SeqIO import parse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class RepresentativeSelector:
    """Class to handle selection of representative ARG sequences"""
    def __init__(self, fasta_file: Path, threads: int = 1):
        self.fasta_file = Path(fasta_file)
        self.threads = threads
        self.logger = logging.getLogger(__name__)
        self._validate_input()
        self.sequences = self._load_sequences()
        
    def _validate_input(self) -> None:
        """Validate input files exist"""
        if not self.fasta_file.exists():
            raise FileNotFoundError(f"FASTA file not found: {self.fasta_file}")

    def _load_sequences(self) -> Dict[str, SeqRecord]:
        """Load and validate protein sequences"""
        sequences = {}
        try:
            for record in parse(str(self.fasta_file), "fasta"):
                # Validate protein sequence
                if not all(c in "ACDEFGHIKLMNPQRSTVWY*X" for c in str(record.seq).upper()):
                    self.logger.warning(f"Sequence {record.id} contains non-protein characters")
                    continue
                sequences[record.id] = record
            
            if not sequences:
                raise ValueError("No valid protein sequences found in input file")
                
            return sequences
        except Exception as e:
            self.logger.error(f"Failed to load sequences: {e}")
            raise

    def _setup_blast_db(self) -> None:
        """Setup BLAST database for input sequences"""
        self.logger.info("Setting up BLAST database...")
        
        # Create a temporary file with validated sequences
        temp_fasta = self.fasta_file.parent / "temp_validated.fasta"

        self.blast_db_prefix = self.fasta_file.parent / "blast_db"

        with open(temp_fasta, 'w') as f:
            for record in self.sequences.values():
                SeqIO.write(record, f, "fasta")
        
        cmd = [
            'makeblastdb',
            '-in', str(temp_fasta),
            '-dbtype', 'prot',
            '-parse_seqids',
            '-out', str(self.blast_db_prefix)
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.debug(f"BLAST DB creation output: {result.stdout}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to create BLAST database: {e.stderr}")
            raise
        finally:
            temp_fasta.unlink(missing_ok=True)

    def _run_blast_comparison(self, output_path: Path, identity: float,
                            coverage: float, evalue: str) -> None:
        """Run BLAST comparison between sequences"""
        self.logger.info("Running BLAST comparison...")
        blast_output = output_path.parent / "blast_results.tsv"
        
        # Create temporary query file
        temp_query = self.fasta_file.parent / "temp_query.fasta"
        with open(temp_query, 'w') as f:
            for record in self.sequences.values():
                SeqIO.write(record, f, "fasta")
        
        cmd = [
            'blastp',
            '-query', str(temp_query),
            '-db', str(self.blast_db_prefix),
            '-out', str(blast_output),
            '-evalue', evalue,
            '-outfmt', '6 qseqid sseqid pident qcovs',
            '-num_threads', str(self.threads),
            '-max_target_seqs', '5000'  # Increase max targets
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.debug(f"BLAST comparison output: {result.stdout}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BLAST comparison failed: {e.stderr}")
            raise
        finally:
            temp_query.unlink(missing_ok=True)

    def _select_representatives(self, output_path: Path) -> None:
        """Select representative sequences from clusters"""
        self.logger.info("Selecting representative sequences...")
        try:
            # Parse BLAST results and create clusters
            clusters = parse_clusters(str(output_path.parent / "blast_results.tsv"))
            
            # Merge overlapping clusters
            merged = merge_overlapping_clusters(clusters)
            
            # Load sequences
            sequences = {record.id: str(record.seq)
                       for record in SeqIO.parse(self.fasta_file, "fasta")}
            
            # Load CARD IDs if available
            card_ids = set()
            card_ids_file = output_path.parent / "card_ids.txt"
            if card_ids_file.exists():
                card_ids = load_card_ids(str(card_ids_file))
            
            # Select representatives
            representatives = select_representatives(merged, sequences, card_ids)
            
            # Write output
            with open(output_path, 'w') as f:
                for rep_id in representatives:
                    f.write(f">{rep_id}\n{sequences[rep_id]}\n")
                    
        except Exception as e:
            self.logger.error(f"Failed to select representatives: {e}")
            raise

    def run_pipeline(self, output_path: Path, identity: float = 90.0, 
                    coverage: float = 70.0, evalue: str = '1e-10') -> None:
        """Run the complete pipeline"""
        self.logger.info("Starting pipeline...")
        
        # Setup blast database
        self._setup_blast_db()
        
        # Run comparisons
        self._run_blast_comparison(output_path, identity, coverage, evalue)
        
        # Select representatives
        self._select_representatives(output_path)
        
        self.logger.info("Pipeline completed successfully")

def validate_files(clusters: Path, fasta: Path, card_ids: Path) -> None:
    """Validate input files exist and are readable"""
    # Check fasta file first as it's the main requirement
    if not fasta.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta}")
    if not fasta.is_file():
        raise ValueError(f"Not a file: {fasta}")
    
    # Create directories for output files if they don't exist
    clusters.parent.mkdir(parents=True, exist_ok=True)
    card_ids.parent.mkdir(parents=True, exist_ok=True)
    
    # For clusters and card_ids, we'll create them if they don't exist
    if not clusters.exists():
        clusters.touch()
    if not card_ids.exists():
        card_ids.touch()

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Select representative ARGs from clusters')
    parser.add_argument('-c', '--clusters', required=True,
                       help='Input cluster file (CD-HIT format)')
    parser.add_argument('-f', '--fasta', required=True,
                       help='Input FASTA file with all sequences')
    parser.add_argument('-o', '--output', required=True,
                       help='Output FASTA file for representative sequences')
    parser.add_argument('--card-ids', required=True,
                       help='File containing CARD database IDs')
    parser.add_argument('--threads', type=int, default=1,
                       help='Number of threads for BLAST')
    parser.add_argument('--identity', type=float, default=90.0,
                       help='Minimum identity threshold')
    parser.add_argument('--coverage', type=float, default=70.0,
                       help='Minimum coverage threshold')
    parser.add_argument('--evalue', type=str, default='1e-10',
                       help='E-value threshold')
    return parser.parse_args()

def load_card_ids(filename: str) -> Set[str]:
    """Load CARD database IDs from file."""
    with open(filename) as f:
        return set(line.strip() for line in f)

def parse_clusters(cluster_file: str) -> Dict[str, Set[str]]:
    """Parse CD-HIT cluster file into groups of sequence IDs."""
    clusters = defaultdict(set)
    current_cluster = None
    
    with open(cluster_file) as f:
        for line in f:
            if line.startswith('>Cluster'):
                current_cluster = line.strip().split()[1]
            elif line.strip():
                seq_id = line.strip().split('>')[1].split('...')[0]
                clusters[current_cluster].add(seq_id)
    
    return clusters

def merge_overlapping_clusters(clusters: Dict[str, Set[str]]) -> List[Set[str]]:
    """Merge clusters that share common sequences."""
    merged = []
    processed = set()
    
    for cluster_id, members in clusters.items():
        if cluster_id in processed:
            continue
            
        current_group = members.copy()
        processed.add(cluster_id)
        
        # Check for overlaps with other clusters
        changed = True
        while changed:
            changed = False
            for other_id, other_members in clusters.items():
                if other_id not in processed and current_group & other_members:
                    current_group |= other_members
                    processed.add(other_id)
                    changed = True
                    
        merged.append(current_group)
    
    return merged

def select_representatives(merged_clusters: List[Set[str]], 
                         sequences: Dict[str, str],
                         card_ids: Set[str]) -> List[str]:
    """Select representative sequence for each cluster, preferring CARD sequences."""
    representatives = []
    
    for cluster in merged_clusters:
        # First try to find a CARD sequence
        card_matches = cluster & card_ids
        if card_matches:
            representatives.append(next(iter(card_matches)))
        else:
            # If no CARD sequence, take the first one
            representatives.append(next(iter(cluster)))
    
    return representatives

def run_blast_comparison(fasta_file: str, output_file: str, threads: int,
                        evalue: str) -> None:
    """Run BLAST all-vs-all comparison"""
    import subprocess
    
    # Create BLAST database
    subprocess.run(['makeblastdb', '-in', fasta_file, '-dbtype', 'prot'], check=True)
    
    # Run BLAST
    cmd = [
        'blastp',
        '-query', fasta_file,
        '-db', fasta_file,
        '-out', output_file,
        '-evalue', evalue,
        '-outfmt', '6 qseqid sseqid pident qcovs',
        '-num_threads', str(threads)
    ]
    subprocess.run(cmd, check=True)

def extract_card_ids(fasta_file: str, output_file: str) -> None:
    """Extract sequence IDs from FASTA file"""
    with open(fasta_file) as f, open(output_file, 'w') as out:
        for line in f:
            if line.startswith('>'):
                out.write(line[1:].strip() + '\n')

def main():
    """Main entry point"""
    args = parse_args()
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    try:
        # Convert to Path objects
        clusters_path = Path(args.clusters)
        fasta_path = Path(args.fasta)
        card_ids_path = Path(args.card_ids)
        output_path = Path(args.output)
        
        # Create output directory if it doesn't exist
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Validate input files
        validate_files(clusters_path, fasta_path, card_ids_path)
        
        # Initialize selector and run pipeline
        selector = RepresentativeSelector(fasta_path, args.threads)
        selector.run_pipeline(
            output_path,
            args.identity,
            args.coverage,
            args.evalue
        )
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        raise

if __name__ == "__main__":
    main()
