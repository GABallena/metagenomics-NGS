from dataclasses import dataclass
from pathlib import Path
import logging
import pysam
import numpy as np
from typing import Tuple, Optional

@dataclass
class Region:
    """Represents a genomic region"""
    chrom: str
    start: int
    end: int
    
    def __iter__(self):
        return iter((self.chrom, self.start, self.end))

class MobilomeAnalyzer:
    """Analyzes mobilome characteristics in genomic regions"""
    def __init__(self, bam_file: Path, **params):
        self.bam_file = Path(bam_file)
        self.params = self._validate_params(params)
        self.logger = self._setup_logging()
        self._validate_bam_file()

    def _validate_params(self, params: dict) -> dict:
        """Validate and set default parameters"""
        defaults = {
            'min_coverage': 0.7,
            'max_soft_clip': 10,
            'max_insert_size': 1050
        }
        return {**defaults, **params}

    def _setup_logging(self):
        """Setup logging"""
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        return logger

    def _validate_bam_file(self):
        """Check if BAM file is valid and indexed"""
        try:
            with pysam.AlignmentFile(self.bam_file, "rb") as bam:
                if not bam.check_index():
                    self.logger.warning("BAM file is not indexed, performance may be impacted")
        except Exception as e:
            raise ValueError(f"Invalid BAM file: {e}")

    def check_conjunction_reads(self, arg_region: Tuple[str, int, int], neighbor_region: Tuple[str, int, int]) -> bool:
        """
        Check for conjunction reads between ARG and neighboring gene
        
        Args:
            arg_region: Region coordinates for ARG
            neighbor_region: Region coordinates for neighboring gene
        
        Returns:
            bool: True if valid conjunction reads found
        """
        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            reads = []
            for read in bam.fetch(arg_region[0], arg_region[1], arg_region[2]):
                if (read.mapping_quality > 0 and 
                    read.get_tag("AS")/len(read.query_sequence) > 0.2 and
                    len(read.get_soft_clips()) < self.params['max_soft_clip']):
                    reads.append(read)
                    
            # Check paired reads
            for read in reads:
                if (read.is_paired and 
                    read.template_length < self.params['max_insert_size'] and
                    self.is_mate_in_region(read, neighbor_region)):
                    return True
        return False
    
    def is_mate_in_region(self, read, region: Tuple[str, int, int]) -> bool:
        """Check if mate is in specified region"""
        return (read.next_reference_start >= region[1] and 
                read.next_reference_start <= region[2])
    
    def calculate_coverage(self, region: Tuple[str, int, int]) -> float:
        """Calculate coverage for a given region"""
        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            coverage = bam.count_coverage(region[0], region[1], region[2])
            avg_coverage = np.mean([sum(x) for x in coverage])
            return avg_coverage / (region[2] - region[1])
    
    def analyze_mobilome(self, arg_region: Tuple[str, int, int], neighbor_region: Tuple[str, int, int]) -> bool:
        """
        Analyze mobilome based on given criteria
        
        Returns:
            bool: True if mobilome criteria are satisfied
        """
        coverage = self.calculate_coverage(arg_region)
        if coverage < self.params['min_coverage']:
            return False
            
        return self.check_conjunction_reads(arg_region, neighbor_region)
