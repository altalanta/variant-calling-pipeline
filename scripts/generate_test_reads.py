#!/usr/bin/env python3

"""
Generate synthetic paired-end FASTQ files for testing
Creates small FASTQ files that will align to the test reference
"""

import gzip
import random
import argparse
from pathlib import Path

# Test reference sequence (from make_ref_bundle.sh)
TEST_REF = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 14 + "ATCGATCGATCGATCGATCG"

def generate_read_pair(ref_seq, read_length=75, insert_size=200):
    """Generate a paired-end read from reference sequence"""
    
    # Random start position (ensuring both reads fit)
    max_start = len(ref_seq) - insert_size - read_length
    if max_start <= 0:
        max_start = len(ref_seq) - read_length
        insert_size = read_length
    
    start_pos = random.randint(0, max_start)
    
    # Extract forward read
    read1_seq = ref_seq[start_pos:start_pos + read_length]
    
    # Extract reverse read (reverse complement)
    read2_start = start_pos + insert_size - read_length
    if read2_start + read_length > len(ref_seq):
        read2_start = len(ref_seq) - read_length
    
    read2_seq = ref_seq[read2_start:read2_start + read_length]
    read2_seq = reverse_complement(read2_seq)
    
    return read1_seq, read2_seq

def reverse_complement(seq):
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def introduce_errors(seq, error_rate=0.01):
    """Introduce random sequencing errors"""
    bases = ['A', 'T', 'C', 'G']
    seq_list = list(seq)
    
    for i in range(len(seq_list)):
        if random.random() < error_rate:
            # Choose different base
            new_bases = [b for b in bases if b != seq_list[i]]
            seq_list[i] = random.choice(new_bases)
    
    return ''.join(seq_list)

def generate_quality_scores(length, min_qual=20, max_qual=40):
    """Generate random quality scores in FASTQ format"""
    # Quality scores: 33 offset, so chr(33+20) = '5', chr(33+40) = 'I'
    return ''.join(chr(33 + random.randint(min_qual, max_qual)) for _ in range(length))

def generate_fastq_files(output_dir, sample_id="NA00001", num_reads=1000, read_length=75):
    """Generate paired FASTQ files"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    read1_file = output_dir / f"tiny_reads_1.fq.gz"
    read2_file = output_dir / f"tiny_reads_2.fq.gz"
    
    with gzip.open(read1_file, 'wt') as f1, gzip.open(read2_file, 'wt') as f2:
        
        for i in range(num_reads):
            # Generate read pair
            read1_seq, read2_seq = generate_read_pair(TEST_REF, read_length)
            
            # Introduce some sequencing errors
            read1_seq = introduce_errors(read1_seq)
            read2_seq = introduce_errors(read2_seq)
            
            # Generate quality scores
            qual1 = generate_quality_scores(len(read1_seq))
            qual2 = generate_quality_scores(len(read2_seq))
            
            # Write FASTQ records
            read_id = f"@{sample_id}_read_{i+1:06d}"
            
            # Read 1
            f1.write(f"{read_id}/1\n")
            f1.write(f"{read1_seq}\n")
            f1.write("+\n")
            f1.write(f"{qual1}\n")
            
            # Read 2
            f2.write(f"{read_id}/2\n")
            f2.write(f"{read2_seq}\n")
            f2.write("+\n")
            f2.write(f"{qual2}\n")
    
    print(f"Generated {num_reads} read pairs:")
    print(f"  - {read1_file}")
    print(f"  - {read2_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate test FASTQ files")
    parser.add_argument("--output-dir", default="tests/smoke", 
                       help="Output directory for FASTQ files")
    parser.add_argument("--sample-id", default="NA00001",
                       help="Sample identifier")
    parser.add_argument("--num-reads", type=int, default=1000,
                       help="Number of read pairs to generate")
    parser.add_argument("--read-length", type=int, default=75,
                       help="Read length")
    
    args = parser.parse_args()
    
    generate_fastq_files(
        args.output_dir,
        args.sample_id,
        args.num_reads,
        args.read_length
    )

if __name__ == "__main__":
    main()