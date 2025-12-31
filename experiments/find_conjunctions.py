import pysam
import argparse

def check_read_alignment(read):
    """Check if read meets basic alignment criteria.

    Criteria (template defaults):
    - aligned_length / read_length > 0.20
    - total soft-clipped bases < 10
    """
    if read.query_length is None or read.query_length == 0:
        return False
    if read.cigartuples is None:
        return False

    aligned_length = sum(length for op, length in read.cigartuples if op == 0)  # M
    alignment_ratio = aligned_length / read.query_length

    soft_clips = sum(length for op, length in read.cigartuples if op == 4)  # S
    return alignment_ratio > 0.20 and soft_clips < 10

def find_conjunction_reads(bam_file, output_file):
    """Find read pairs that connect ARGs to neighboring genes"""
    bam = pysam.AlignmentFile(bam_file, "rb")
    # Note: BAM should be indexed for efficient fetch; this script iterates sequentially.
    
    valid_pairs = set()
    for read in bam.fetch():
        if not read.is_paired or read.is_secondary:
            continue
            
        # Check if read pair meets criteria
        if (check_read_alignment(read) and 
            abs(read.template_length) < 1050):  # 350 × 3 bp
            valid_pairs.add(read.query_name)
    
    # Write results
    with open(output_file, 'w') as out:
        for read_name in valid_pairs:
            out.write(f"{read_name}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    
    find_conjunction_reads(args.bam, args.output)
