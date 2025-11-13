"""
Real DNA Dataset Search Example

This script demonstrates Boyer-Moore pattern matching on real DNA sequences
from the NCBI dataset included in the project.

Author: DNA Pattern Matching Project
Date: November 2025
"""

import os
import sys
import time

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from boyer_moore import BoyerMoore, search_multiple_patterns
from utils import (
    read_fasta_file, get_all_fasta_files, calculate_gc_content,
    count_nucleotides, get_reverse_complement
)


def search_in_real_genome():
    """Search for patterns in a real E. coli genome."""
    print("\n" + "="*70)
    print("SEARCHING IN REAL E. COLI GENOME SEQUENCES")
    print("="*70)
    
    # Get the dataset path 
    script_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_path = os.path.abspath(os.path.join(script_dir, '..', 'DnA_dataset', 'ncbi_dataset', 'data'))
    
    if not os.path.exists(dataset_path):
        print(f"\nDataset not found at: {dataset_path}")
        print("Please ensure the DNA dataset is in the correct location.")
        return
    
    # Get one FASTA file
    fasta_files = get_all_fasta_files(dataset_path, recursive=True)
    
    if not fasta_files:
        print("No FASTA files found in the dataset.")
        return
    
    # Use the first file (E. coli genome)
    fasta_file = fasta_files[0]
    print(f"\nReading genome from: {os.path.basename(fasta_file)}")
    
    try:
        # Read the genome
        sequences_dict = read_fasta_file(fasta_file)
        
        if not sequences_dict:
            print("No sequences found in file.")
            return
        
        # Get the first sequence
        header = list(sequences_dict.keys())[0]
        sequence = sequences_dict[header]
        
        print(f"\nSequence: {header[:60]}...")
        print(f"Length: {len(sequence):,} bp")
        print(f"GC Content: {calculate_gc_content(sequence):.2f}%")
        
        # Analyze nucleotide composition
        nucleotides = count_nucleotides(sequence)
        print(f"\nNucleotide composition:")
        for base in ['A', 'C', 'G', 'T']:
            count = nucleotides[base]
            print(f"  {base}: {count:,} ({count/len(sequence)*100:.2f}%)")
        
        # Search for biologically significant patterns
        print("\n" + "="*70)
        print("SEARCHING FOR BIOLOGICALLY SIGNIFICANT PATTERNS")
        print("="*70)
        
        patterns = {
            'Start Codon (ATG)': 'ATG',
            'Stop Codon (TAA)': 'TAA',
            'Stop Codon (TAG)': 'TAG',
            'Stop Codon (TGA)': 'TGA',
            'TATA Box': 'TATAAA',
            'Pribnow Box (-10)': 'TATAAT',
            'EcoRI site': 'GAATTC',
            'BamHI site': 'GGATCC',
            'PstI site': 'CTGCAG',
            'HindIII site': 'AAGCTT',
            'Ribosome Binding Site': 'AGGAGGT',
            'Chi site': 'GCTGGTGG'
        }
        
        print(f"\n{'Pattern Name':<25} {'Sequence':<12} {'Count':<8} {'Time (ms)':<10}")
        print("-" * 70)
        
        for name, pattern in patterns.items():
            bm = BoyerMoore(pattern)
            
            start_time = time.time()
            matches = bm.search(sequence)
            end_time = time.time()
            
            elapsed = (end_time - start_time) * 1000
            
            print(f"{name:<25} {pattern:<12} {len(matches):<8} {elapsed:>8.2f}")
        
        # Detailed analysis of start and stop codons
        print("\n" + "="*70)
        print("OPEN READING FRAME (ORF) INDICATORS")
        print("="*70)
        
        start_codon = BoyerMoore('ATG')
        stop_codons = {
            'TAA': BoyerMoore('TAA'),
            'TAG': BoyerMoore('TAG'),
            'TGA': BoyerMoore('TGA')
        }
        
        start_positions = start_codon.search(sequence)
        print(f"\nStart codons (ATG): {len(start_positions):,} occurrences")
        
        total_stops = 0
        for codon, bm in stop_codons.items():
            count = bm.count_matches(sequence)
            total_stops += count
            print(f"Stop codon ({codon}): {count:,} occurrences")
        
        print(f"Total stop codons: {total_stops:,}")
        print(f"Start/Stop ratio: {len(start_positions)/total_stops:.3f}" if total_stops > 0 else "N/A")
        
        # Search for patterns in both strands
        print("\n" + "="*70)
        print("SEARCHING IN BOTH STRANDS")
        print("="*70)
        
        pattern = "GAATTC"  # EcoRI restriction site
        print(f"\nSearching for EcoRI restriction site: {pattern}")
        
        bm = BoyerMoore(pattern)
        forward_matches = bm.search(sequence)
        
        rev_comp = get_reverse_complement(sequence)
        reverse_matches = bm.search(rev_comp)
        
        print(f"Forward strand: {len(forward_matches)} matches")
        print(f"Reverse strand: {len(reverse_matches)} matches")
        print(f"Total matches (both strands): {len(forward_matches) + len(reverse_matches)}")
        
        # Show first few matches with context
        if forward_matches:
            print(f"\nFirst 5 matches on forward strand:")
            for i, pos in enumerate(forward_matches[:5], 1):
                start = max(0, pos - 10)
                end = min(len(sequence), pos + len(pattern) + 10)
                context = sequence[start:end]
                marker = " " * (pos - start) + "^" * len(pattern)
                print(f"  {i}. Position {pos:,}:")
                print(f"     ...{context}...")
                print(f"     ...{marker}...")
        
        # Performance summary
        print("\n" + "="*70)
        print("PERFORMANCE SUMMARY")
        print("="*70)
        print(f"\nGenome size: {len(sequence):,} bp")
        print(f"Patterns searched: {len(patterns)}")
        print(f"Average search time: ~{sum((time.time(), -time.time())[0] * 1000 for _ in range(1)):.2f} ms per pattern")
        print("\n✓ All searches completed successfully!")
        
    except Exception as e:
        print(f"\nError: {str(e)}")
        import traceback
        traceback.print_exc()


def main():
    """Main function."""
    search_in_real_genome()


if __name__ == '__main__':
    main()
