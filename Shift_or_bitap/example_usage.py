"""
Example Usage of Shift-Or/Bitap Algorithm Implementation

This script demonstrates how to use all three variants of the Shift-Or algorithm:
1. Exact matching (≤64 bp)
2. Approximate matching (≤64 bp, k=1,2,3)
3. Extended matching (>64 bp)

Author: DNA Pattern Matching Project
Date: November 2025
"""

from shift_or_exact import ShiftOrExact, search_multiple_patterns as search_exact_multi
from shift_or_approximate import ShiftOrApproximate, search_multiple_patterns as search_approx_multi
from shift_or_extended import ShiftOrExtended, search_multiple_patterns as search_extended_multi
from shift_or_utils import read_fasta_file, generate_random_dna


def example_1_exact_matching():
    """Example 1: Basic exact matching."""
    print("=" * 70)
    print("EXAMPLE 1: Exact Matching (≤64 bp)")
    print("=" * 70)

    # Create sample DNA sequence
    text = "ACGTACGTTAGCTAGCTAGCTACGTACGTACGT"
    pattern = "ACGT"

    # Initialize matcher
    matcher = ShiftOrExact(pattern)

    # Search for pattern
    matches = matcher.search(text)

    print(f"Text:    {text}")
    print(f"Pattern: {pattern}")
    print(f"Matches found at positions: {matches}")
    print(f"Number of matches: {len(matches)}\n")

    # Show detailed metrics
    metrics = matcher.search_with_metrics(text)
    print("Detailed Metrics:")
    print(f"  Preprocessing time: {matcher.get_preprocessing_time():.4f} ms")
    print(f"  Search time: {metrics['search_time_ms']:.4f} ms")
    print(f"  Bit operations: {metrics['bit_operations']}")
    print(f"  State vectors: {metrics['state_vectors']}")
    print()


def example_2_approximate_matching():
    """Example 2: Approximate matching with errors."""
    print("=" * 70)
    print("EXAMPLE 2: Approximate Matching (k errors allowed)")
    print("=" * 70)

    text = "ACGTACTTACATAAAGACGTNNNN"
    pattern = "ACGT"

    # Try different error levels
    for k in [1, 2, 3]:
        print(f"\n--- k = {k} (allowing up to {k} error(s)) ---")
        matcher = ShiftOrApproximate(pattern, k=k)
        matches = matcher.search(text)

        print(f"Pattern: {pattern}")
        print(f"Found {len(matches)} matches:")
        for pos, error_level in matches[:5]:  # Show first 5
            match_text = text[pos:pos+len(pattern)]
            print(f"  Position {pos}: '{match_text}' (error level: {error_level})")


def example_3_extended_matching():
    """Example 3: Extended matching for long patterns."""
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Extended Matching (>64 bp patterns)")
    print("=" * 70)

    # Create a long pattern (100 bp)
    pattern_100 = "ACGT" * 25  # 100 bp
    text = "N" * 500 + pattern_100 + "N" * 500

    print(f"Pattern length: {len(pattern_100)} bp")
    print(f"Text length: {len(text)} bp")

    matcher = ShiftOrExtended(pattern_100)
    matches = matcher.search(text)

    print(f"Matches found at positions: {matches}")
    print(f"Number of matches: {len(matches)}")

    # Show metrics
    metrics = matcher.search_with_metrics(text)
    print(f"\nMulti-word implementation metrics:")
    print(f"  State vectors (words): {metrics['state_vectors']}")
    print(f"  Bit operations: {metrics['bit_operations']}")
    print()


def example_4_multiple_patterns():
    """Example 4: Searching for multiple patterns."""
    print("=" * 70)
    print("EXAMPLE 4: Multiple Pattern Search")
    print("=" * 70)

    text = generate_random_dna(1000, seed=42)
    patterns = ["ACGT", "TTTT", "GGGG", "CCCC", "ATCG"]

    print(f"Searching for {len(patterns)} patterns in {len(text)} bp sequence...")

    # Exact matching
    print("\n--- Exact Matching ---")
    results_exact = search_exact_multi(text, patterns)
    for pattern, matches in results_exact.items():
        print(f"  {pattern}: {len(matches)} matches")

    # Approximate matching (k=1)
    print("\n--- Approximate Matching (k=1) ---")
    results_approx = search_approx_multi(text, patterns, k=1)
    for pattern, matches in results_approx.items():
        print(f"  {pattern}: {len(matches)} matches")


def example_5_real_fasta():
    """Example 5: Working with real FASTA files."""
    print("\n" + "=" * 70)
    print("EXAMPLE 5: Processing FASTA Files")
    print("=" * 70)

    print("\nTo use with your own FASTA file:")
    print("""
# Read FASTA file
from shift_or_utils import read_fasta_file

sequences = read_fasta_file("your_genome.fna")
text = list(sequences.values())[0]  # Get first sequence

# Search for pattern
pattern = "ACGTACGT"
matcher = ShiftOrExact(pattern)
matches = matcher.search(text)

print(f"Found {len(matches)} matches in genome")
""")

    # Demonstrate with generated sequence
    print("\nDemonstration with generated 10,000 bp sequence:")
    text = generate_random_dna(10000, seed=42)
    pattern = "ACGTACGT"

    matcher = ShiftOrExact(pattern)
    matches = matcher.search(text)

    print(f"Pattern: {pattern}")
    print(f"Sequence length: {len(text):,} bp")
    print(f"Matches found: {len(matches)}")
    print(f"First 10 match positions: {matches[:10]}")


def example_6_performance_comparison():
    """Example 6: Performance comparison across variants."""
    print("\n" + "=" * 70)
    print("EXAMPLE 6: Performance Comparison")
    print("=" * 70)

    import time

    # Generate test data
    text = generate_random_dna(100000, seed=42)
    pattern_short = "ACGTACGT"  # 8 bp
    pattern_64 = "ACGT" * 16     # 64 bp
    pattern_100 = "ACGT" * 25    # 100 bp

    print(f"Testing on {len(text):,} bp sequence\n")

    # Exact matching - short pattern
    start = time.perf_counter()
    matcher1 = ShiftOrExact(pattern_short)
    matches1 = matcher1.search(text)
    time1 = (time.perf_counter() - start) * 1000
    print(f"Exact (8 bp):    {time1:.3f} ms, {len(matches1)} matches")

    # Exact matching - 64 bp pattern
    start = time.perf_counter()
    matcher2 = ShiftOrExact(pattern_64)
    matches2 = matcher2.search(text)
    time2 = (time.perf_counter() - start) * 1000
    print(f"Exact (64 bp):   {time2:.3f} ms, {len(matches2)} matches")

    # Extended matching - 100 bp pattern
    start = time.perf_counter()
    matcher3 = ShiftOrExtended(pattern_100)
    matches3 = matcher3.search(text)
    time3 = (time.perf_counter() - start) * 1000
    print(f"Extended (100 bp): {time3:.3f} ms, {len(matches3)} matches")

    # Approximate matching - k=1
    start = time.perf_counter()
    matcher4 = ShiftOrApproximate(pattern_short, k=1)
    matches4 = matcher4.search(text)
    time4 = (time.perf_counter() - start) * 1000
    print(f"Approximate k=1 (8 bp): {time4:.3f} ms, {len(matches4)} matches")

    # Approximate matching - k=3
    start = time.perf_counter()
    matcher5 = ShiftOrApproximate(pattern_short, k=3)
    matches5 = matcher5.search(text)
    time5 = (time.perf_counter() - start) * 1000
    print(f"Approximate k=3 (8 bp): {time5:.3f} ms, {len(matches5)} matches")


def main():
    """Run all examples."""
    print("\n" + "=" * 70)
    print("SHIFT-OR/BITAP ALGORITHM - USAGE EXAMPLES")
    print("=" * 70 + "\n")

    example_1_exact_matching()
    example_2_approximate_matching()
    example_3_extended_matching()
    example_4_multiple_patterns()
    example_5_real_fasta()
    example_6_performance_comparison()

    print("\n" + "=" * 70)
    print("Examples completed! Check the code for details.")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
