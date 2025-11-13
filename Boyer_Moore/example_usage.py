"""
Example Usage of Boyer-Moore Algorithm for DNA Pattern Matching

This script demonstrates various use cases and features of the Boyer-Moore
implementation for DNA sequence analysis.
"""

import os
import sys

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from boyer_moore import BoyerMoore, search_multiple_patterns, find_approximate_matches
from utils import (
    read_fasta_file, generate_random_dna, get_reverse_complement,
    calculate_gc_content, validate_dna_sequence, count_nucleotides
)


def example_1_basic_search():
    """Example 1: Basic pattern search in a DNA sequence."""
    print("\n" + "="*70)
    print("EXAMPLE 1: Basic Pattern Search")
    print("="*70)
    
    # Create a sample DNA sequence
    text = "ACGTACGTATGCATGCACGTATGCTAGC"
    pattern = "ATGC"
    
    print(f"Text:    {text}")
    print(f"Pattern: {pattern}")
    
    # Create Boyer-Moore instance and search
    bm = BoyerMoore(pattern)
    matches = bm.search(text)
    
    print(f"\nPattern found at positions: {matches}")
    print(f"Total matches: {len(matches)}")
    
    # Show matches in context
    for pos in matches:
        start = max(0, pos - 5)
        end = min(len(text), pos + len(pattern) + 5)
        context = text[start:end]
        marker = " " * (pos - start) + "^" * len(pattern)
        print(f"\nPosition {pos}: ...{context}...")
        print(f"             ...{marker}...")


def example_2_restriction_enzymes():
    """Example 2: Finding restriction enzyme sites."""
    print("\n" + "="*70)
    print("EXAMPLE 2: Finding Restriction Enzyme Sites")
    print("="*70)
    
    # Generate a random DNA sequence
    sequence = generate_random_dna(1000, seed=42)
    
    # Common restriction enzyme recognition sites
    restriction_sites = {
        'EcoRI': 'GAATTC',
        'BamHI': 'GGATCC',
        'PstI': 'CTGCAG',
        'HindIII': 'AAGCTT',
        'SalI': 'GTCGAC',
        'XhoI': 'CTCGAG'
    }
    
    print(f"Searching in a {len(sequence)} bp sequence...")
    print()
    
    for enzyme_name, site in restriction_sites.items():
        bm = BoyerMoore(site)
        matches = bm.search(sequence)
        
        if matches:
            print(f"✓ {enzyme_name:10s} ({site}): {len(matches):3d} sites at positions {matches[:5]}"
                  + ("..." if len(matches) > 5 else ""))
        else:
            print(f"✗ {enzyme_name:10s} ({site}): No sites found")


def example_3_multiple_patterns():
    """Example 3: Searching for multiple patterns simultaneously."""
    print("\n" + "="*70)
    print("EXAMPLE 3: Multiple Pattern Search")
    print("="*70)
    
    text = "ACGTACGTATGCATGCACGTATGCTAGCATGCACGT"
    patterns = ["ATGC", "ACGT", "CGTA", "TAGC"]
    
    print(f"Text: {text}")
    print(f"Searching for {len(patterns)} patterns: {patterns}")
    print()
    
    results = search_multiple_patterns(text, patterns)
    
    for pattern, positions in results.items():
        print(f"Pattern '{pattern}': {len(positions)} matches at {positions}")


def example_4_first_occurrence():
    """Example 4: Finding only the first occurrence."""
    print("\n" + "="*70)
    print("EXAMPLE 4: Finding First Occurrence")
    print("="*70)
    
    text = "ACGTACGTATGCATGCACGTATGCTAGC"
    pattern = "ATGC"
    
    bm = BoyerMoore(pattern)
    first_pos = bm.search_first(text)
    
    print(f"Text:    {text}")
    print(f"Pattern: {pattern}")
    print(f"\nFirst occurrence at position: {first_pos}")


def example_5_reverse_complement():
    """Example 5: Searching in reverse complement."""
    print("\n" + "="*70)
    print("EXAMPLE 5: Searching in Reverse Complement")
    print("="*70)
    
    sequence = "ACGTACGTATGCATGC"
    pattern = "ATGC"
    
    print(f"Original sequence: {sequence}")
    print(f"Pattern to find:   {pattern}")
    
    # Search in forward strand
    bm = BoyerMoore(pattern)
    forward_matches = bm.search(sequence)
    
    # Search in reverse complement
    rev_comp = get_reverse_complement(sequence)
    reverse_matches = bm.search(rev_comp)
    
    print(f"\nReverse complement: {rev_comp}")
    print(f"\nForward strand matches: {forward_matches} ({len(forward_matches)} total)")
    print(f"Reverse strand matches: {reverse_matches} ({len(reverse_matches)} total)")


def example_6_sequence_analysis():
    """Example 6: Sequence analysis with pattern matching."""
    print("\n" + "="*70)
    print("EXAMPLE 6: Sequence Analysis with Pattern Matching")
    print("="*70)
    
    # Generate a sequence
    sequence = generate_random_dna(500, seed=123)
    
    # Analyze the sequence
    print(f"Sequence length: {len(sequence)} bp")
    print(f"GC content: {calculate_gc_content(sequence):.2f}%")
    
    nucleotide_counts = count_nucleotides(sequence)
    print(f"\nNucleotide composition:")
    for base, count in nucleotide_counts.items():
        if count > 0:
            print(f"  {base}: {count} ({count/len(sequence)*100:.1f}%)")
    
    # Search for common motifs
    print(f"\nSearching for common motifs:")
    motifs = {
        'CpG island': 'CG',
        'Poly-A tail': 'AAAA',
        'Kozak sequence': 'CCACCATGG'
    }
    
    for motif_name, motif_seq in motifs.items():
        bm = BoyerMoore(motif_seq)
        count = bm.count_matches(sequence)
        print(f"  {motif_name} ({motif_seq}): {count} occurrences")


def example_7_performance_comparison():
    """Example 7: Performance comparison with different pattern lengths."""
    print("\n" + "="*70)
    print("EXAMPLE 7: Performance Comparison")
    print("="*70)
    
    import time
    
    # Generate a large text
    text = generate_random_dna(100000, seed=42)
    print(f"Text length: {len(text):,} bp")
    print()
    
    # Test with different pattern lengths
    pattern_lengths = [5, 10, 20, 50, 100]
    
    print(f"{'Pattern Length':<15} {'Matches':<10} {'Time (ms)':<12} {'Speed (MB/s)':<12}")
    print("-" * 60)
    
    for length in pattern_lengths:
        # Extract pattern from text to ensure matches exist
        pattern = text[1000:1000+length]
        
        # Time the search
        bm = BoyerMoore(pattern)
        start_time = time.time()
        matches = bm.search(text)
        end_time = time.time()
        
        elapsed = end_time - start_time
        speed = len(text) / elapsed / 1_000_000 if elapsed > 0 else 0
        
        print(f"{length:<15} {len(matches):<10} {elapsed*1000:>10.3f}  {speed:>10.2f}")


def example_8_validate_sequence():
    """Example 8: Validating DNA sequences before search."""
    print("\n" + "="*70)
    print("EXAMPLE 8: DNA Sequence Validation")
    print("="*70)
    
    test_sequences = [
        ("ACGTACGT", "Valid DNA sequence"),
        ("ACGTACGTN", "Valid with ambiguous base (N)"),
        ("ACGTXYZ", "Invalid - contains non-DNA characters"),
        ("acgtacgt", "Valid - lowercase (will be normalized)"),
    ]
    
    for seq, description in test_sequences:
        is_valid = validate_dna_sequence(seq)
        status = "✓ VALID" if is_valid else "✗ INVALID"
        print(f"{status:10s} - {description:40s} [{seq}]")


def example_9_statistics():
    """Example 9: Getting algorithm statistics."""
    print("\n" + "="*70)
    print("EXAMPLE 9: Algorithm Statistics and Internals")
    print("="*70)
    
    pattern = "ATGCATGC"
    bm = BoyerMoore(pattern)
    
    stats = bm.get_statistics()
    
    print(f"Pattern: {stats['pattern']}")
    print(f"Pattern Length: {stats['pattern_length']}")
    print(f"Alphabet Size: {stats['alphabet_size']}")
    print(f"Alphabet: {stats['bad_char_table_keys']}")
    print(f"\nGood Suffix Table: {stats['good_suffix_table']}")


def example_10_overlapping_matches():
    """Example 10: Finding overlapping matches."""
    print("\n" + "="*70)
    print("EXAMPLE 10: Overlapping Pattern Matches")
    print("="*70)
    
    # Create a sequence with overlapping patterns
    text = "AAAAAAAAA"
    pattern = "AAA"
    
    print(f"Text:    {text}")
    print(f"Pattern: {pattern}")
    
    bm = BoyerMoore(pattern)
    matches = bm.search(text)
    
    print(f"\nOverlapping matches found at: {matches}")
    print(f"Total matches: {len(matches)}")
    
    # Visualize matches
    print("\nVisualization:")
    for i, pos in enumerate(matches, 1):
        marker = " " * pos + "^" * len(pattern) + f" Match {i}"
        print(f"{text}")
        print(f"{marker}")


def main():
    """Run all examples."""
    print("\n")
    print("="*70)
    print("BOYER-MOORE ALGORITHM - EXAMPLE USAGE DEMONSTRATIONS")
    print("="*70)
    
    examples = [
        example_1_basic_search,
        example_2_restriction_enzymes,
        example_3_multiple_patterns,
        example_4_first_occurrence,
        example_5_reverse_complement,
        example_6_sequence_analysis,
        example_7_performance_comparison,
        example_8_validate_sequence,
        example_9_statistics,
        example_10_overlapping_matches
    ]
    
    for example_func in examples:
        try:
            example_func()
        except Exception as e:
            print(f"\nError in {example_func.__name__}: {str(e)}")
    
    print("\n" + "="*70)
    print("ALL EXAMPLES COMPLETED")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
