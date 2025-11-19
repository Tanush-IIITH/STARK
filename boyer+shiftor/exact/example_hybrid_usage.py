"""
Example Usage of Hybrid DNA Matcher (EXACT MATCHING)

Demonstrates Architecture A with both algorithms doing exact matching:
- Boyer-Moore: Fast with skipping
- Shift-Or: Bit-parallel exact matching

Author: STARK_5 Analysis Group
"""

from hybrid_dna_matcher import HybridDNAMatcher, search_dna_hybrid


def example_1_basic_exact():
    """Basic exact matching example."""
    print("="*70)
    print("EXAMPLE 1: Basic Hybrid Exact Matching")
    print("="*70)

    text = "ACGTACGTACGTAAGTACGTACGTACGTACGT"
    pattern = "ACGTACGT"

    matches, stats = search_dna_hybrid(text, pattern, verbose=True)

    print(f"\nFound {len(matches)} exact matches:")
    for pos in matches:
        print(f"  Position {pos}")


def example_2_high_similarity():
    """Example with high similarity regions to trigger Shift-Or."""
    print("\n\n" + "="*70)
    print("EXAMPLE 2: High Similarity Regions")
    print("="*70)

    # Text with many similar but not exact matches
    text = "ACGTACGAACGTACGTACGTACGGACGTACGT"
    pattern = "ACGTACGT"

    print(f"Text:    {text}")
    print(f"Pattern: {pattern}")
    print()

    matches, stats = search_dna_hybrid(
        text, pattern,
        pmd_threshold=0.75,
        verbose=True
    )

    print(f"\nExact matches: {len(matches)}")


def example_3_threshold_comparison():
    """Compare different PMD thresholds."""
    print("\n\n" + "="*70)
    print("EXAMPLE 3: PMD Threshold Effect")
    print("="*70)

    text = "ACGTACGTACGTACGTAAGTACGTACGTACGT" * 10
    pattern = "ACGTACGT"

    thresholds = [0.60, 0.75, 0.85, 0.95]

    print(f"Text length: {len(text):,} bp")
    print(f"Pattern: {pattern}\n")

    for threshold in thresholds:
        print(f"{'='*50}")
        print(f"PMD Threshold: {threshold}")
        print(f"{'='*50}")

        matches, stats = search_dna_hybrid(
            text, pattern,
            pmd_threshold=threshold,
            verbose=False
        )

        print(f"Matches: {len(matches)}")
        print(f"Boyer-Moore matches: {stats['boyer_moore_matches']}")
        print(f"Shift-Or matches: {stats['shift_or_matches']}")
        print(f"Shift-Or triggers: {stats['shift_or_triggers']}")
        print(f"Positions skipped: {stats['total_skips']:,}")

        if stats['boyer_moore_scans'] > 0:
            trigger_rate = (stats['shift_or_triggers'] / 
                           stats['boyer_moore_scans'] * 100)
            skip_rate = (stats['total_skips'] / 
                        stats['total_positions_scanned'] * 100)
            print(f"Trigger rate: {trigger_rate:.2f}%")
            print(f"Skip rate: {skip_rate:.2f}%")
        print()


def example_4_pattern_lengths():
    """Test with different pattern lengths."""
    print("\n\n" + "="*70)
    print("EXAMPLE 4: Different Pattern Lengths")
    print("="*70)

    import random
    random.seed(42)

    # Generate random DNA sequence
    text = ''.join(random.choice('ACGT') for _ in range(1000))

    pattern_lengths = [8, 16, 32, 64]

    for length in pattern_lengths:
        pattern = ''.join(random.choice('ACGT') for _ in range(length))

        print(f"\nPattern length: {length} bp")
        matches, stats = search_dna_hybrid(text, pattern, verbose=False)

        print(f"Matches: {len(matches)}")
        print(f"Shift-Or triggers: {stats['shift_or_triggers']}")
        print(f"Skip rate: {stats['total_skips']/stats['total_positions_scanned']*100:.1f}%")


def example_5_fasta_integration():
    """Example with FASTA file reading."""
    print("\n\n" + "="*70)
    print("EXAMPLE 5: FASTA File Integration")
    print("="*70)

    def read_fasta(filename):
        """Read DNA sequence from FASTA file."""
        sequence = []
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence.append(line.strip().upper())
        return ''.join(sequence)

    print("Example function to read FASTA files:")
    print("  text = read_fasta('path/to/file.fna')")
    print("  matches, stats = search_dna_hybrid(text, pattern, verbose=True)")
    print("\nUncomment the code in this function to test with real files.")

    # Uncomment to test with actual FASTA file:
    # try:
    #     text = read_fasta('synthetic_datasets/SYNTH_100K/SYNTH_100K.fna')
    #     pattern = "ATCGATCGATCG"
    #     matches, stats = search_dna_hybrid(text, pattern, verbose=True)
    #     print(f"\nFound {len(matches)} matches in FASTA file")
    # except FileNotFoundError:
    #     print("FASTA file not found")


if __name__ == "__main__":
    example_1_basic_exact()
    example_2_high_similarity()
    example_3_threshold_comparison()
    example_4_pattern_lengths()
    example_5_fasta_integration()

    print("\n" + "="*70)
    print("All examples completed!")
    print("="*70)
