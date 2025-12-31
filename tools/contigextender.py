import os
import sys
import logging
from pathlib import Path
from typing import Tuple, List
from dataclasses import dataclass
from tqdm import tqdm
import subprocess
import time
from threading import Thread


@dataclass
class ScaffoldConfig:
    """Configuration for scaffold processing"""
    coverage: int = 50
    min_overlap: int = 30
    extend_tolerance: float = 2.5
    stop_length: int = 250
    threads: int = 8
    complex_threshold: int = -1
    long_reads_dir: str = "long_reads"
    filtered_contigs_dir: str = "filtered_contigs"
    scaffolds_dir: str = "scaffolds"
    logs_dir: str = "logs"


class ScaffoldProcessor:
    """Processor for running scaffold-related tasks"""
    def __init__(self, config: ScaffoldConfig):
        self.config = config
        self.logger = self._setup_logging()

    def _setup_logging(self) -> logging.Logger:
        """Setup logging configuration"""
        os.makedirs(self.config.logs_dir, exist_ok=True)
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
            handlers=[
                logging.FileHandler(f"{self.config.logs_dir}/scaffolding.log"),
                logging.StreamHandler(sys.stdout),
            ],
        )
        return logging.getLogger(__name__)

    def run_shell_command(self, command: str) -> None:
        """Execute shell command with error handling"""
        self.logger.debug(f"Executing command: {command}")
        result = os.system(command)
        if result != 0:
            self.logger.error(f"Command failed with exit code {result}: {command}")
            raise RuntimeError(f"Command failed: {command}")

    def run_contigextender(self, sample: str) -> Tuple[str, str]:
        """Run contig extender for a given sample"""
        self.logger.info(f"Running ContigExtender for sample {sample}")
        try:
            outdir = Path(self.config.scaffolds_dir) / sample
            scaffolds = Path(self.config.filtered_contigs_dir) / f"{sample}_final_filtered_contigs.fa"
            longreads = Path(self.config.long_reads_dir) / f"{sample}_longreads.fastq"
            filled_scaffolds = outdir / "contigextender_scaffolds.fasta"
            gap_stats = outdir / "extension_stats.txt"
            log_file = Path(self.config.logs_dir) / f"contigextender_{sample}.log"

            outdir.mkdir(parents=True, exist_ok=True)

            command = (
                f"contig-extender/dist/extender_wrapper "
                f"--coverage {self.config.coverage} "
                f"--min-overlap-length {self.config.min_overlap} "
                f"--extend-tolerance {self.config.extend_tolerance} "
                f"--stop-length {self.config.stop_length} "
                f"--threads {self.config.threads} "
                f"--complex-threshold {self.config.complex_threshold} "
                f"--out {outdir} {scaffolds} {longreads} 2> {log_file}"
            )
            self.run_shell_command(command)

            if not filled_scaffolds.exists() or not filled_scaffolds.stat().st_size:
                raise FileNotFoundError(f"ContigExtender failed to create output for {sample}")

            stats_command = f"seqkit stats {filled_scaffolds} > {gap_stats}"
            self.run_shell_command(stats_command)

            return str(filled_scaffolds), str(gap_stats)

        except Exception as e:
            self.logger.error(f"Error processing sample {sample}: {e}")
            raise


def refresh_sudo_periodically(interval: int = 600):
    """Keep sudo active by refreshing the timestamp every `interval` seconds."""
    def refresh_loop():
        while True:
            subprocess.run(['sudo', '-v'], check=False)
            time.sleep(interval)

    Thread(target=refresh_loop, daemon=True).start()


def get_samples(filtered_contigs_dir: str) -> List[str]:
    """Discover samples based on the filtered contigs directory"""
    sample_dir = Path(filtered_contigs_dir)
    samples = [
        file.stem.replace("_final_filtered_contigs", "")
        for file in sample_dir.glob("*_final_filtered_contigs.fa")
    ]
    if not samples:
        raise FileNotFoundError(f"No samples found in the {filtered_contigs_dir} directory.")
    return samples


def main():
    """Main function to process scaffolds"""
    # Ensure sudo privileges upfront
    try:
        subprocess.run(['sudo', '-v'], check=True)
    except subprocess.CalledProcessError:
        print("Failed to obtain sudo privileges. Exiting.")
        sys.exit(1)

    # Start refreshing sudo timestamp in the background
    refresh_sudo_periodically(interval=600)

    config = ScaffoldConfig()
    processor = ScaffoldProcessor(config)

    try:
        samples = get_samples(config.filtered_contigs_dir)
        processor.logger.info(f"Found {len(samples)} samples to process.")

        for sample in tqdm(samples, desc="Processing samples"):
            processor.logger.info(f"Starting processing for sample {sample}")
            try:
                scaffolds, stats = processor.run_contigextender(sample)
                processor.logger.info(f"Completed processing for sample {sample}")
                processor.logger.debug(f"Scaffolds: {scaffolds}, Stats: {stats}")
            except Exception as e:
                processor.logger.error(f"Failed processing for sample {sample}: {e}")
                continue

    except Exception as e:
        processor.logger.error(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
