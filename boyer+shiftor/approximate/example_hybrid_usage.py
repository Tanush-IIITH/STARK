"""
Example Usage of Hybrid DNA Matcher (Updated)

Demonstrates Architecture A using existing Boyer-Moore class
"""

from hybrid_dna_matcher import HybridDNAMatcher, search_dna_hybrid


def example_1_basic_usage():
    """Basic usage with your Boyer-Moore class."""
    print("="*70)
    print("EXAMPLE 1: Basic Hybrid Search")
    print("="*70)

    text = "ACGTACGTACGTAAGTACGTACGTACGTACGT"
    pattern = "ACGTACGT"

    matches, stats = search_dna_hybrid(text, pattern, verbose=True)

    print(f"\nFound {len(matches)} matches:")
    for pos, match_type, errors in matches:
        print(f"  Position {pos}: {match_type} (errors={errors})")


def example_2_with_mutations():
    """Example with mutations to test Shift-Or trigger."""
    print("\n\n" + "="*70)
    print("EXAMPLE 2: Search with Mutations (PMD Trigger)")
    print("="*70)

    # Text with 1-error mutations
    text = "ACGTACGGGACGTACGTAAGTACGTCCGTACGT"  # Mutations at pos 7, 26
    pattern = "ACGTACGT"

    matches, stats = search_dna_hybrid(
        text, pattern, 
        pmd_threshold=0.75, 
        k_errors=1, 
        verbose=True
    )

    print(f"\nMatches: {len(matches)}")


def example_3_threshold_comparison():
    """Compare different PMD thresholds."""
    print("\n\n" + "="*70)
    print("EXAMPLE 3: PMD Threshold Comparison")
    print("="*70)

    text = "ACGTACGGGACGTACGTAAGTACGTCCGTACGT" * 5
    pattern = "ACGTACGT"

    thresholds = [0.60, 0.75, 0.85, 0.95]

    for threshold in thresholds:
        print(f"\n--- PMD Threshold: {threshold} ---")
        matches, stats = search_dna_hybrid(
            text, pattern,
            pmd_threshold=threshold,
            k_errors=1,
            verbose=False
        )

        print(f"Matches: {len(matches)} "
              f"(Exact: {stats['exact_matches']}, "
              f"Approx: {stats['approximate_matches']})")
        print(f"Shift-Or triggers: {stats['shift_or_triggers']}")
        print(f"Trigger rate: {stats['shift_or_triggers']/stats['boyer_moore_scans']*100:.1f}%")


def example_4_read_fasta():
    """Example reading from FASTA file."""
    print("\n\n" + "="*70)
    print("EXAMPLE 4: FASTA File Search")
    print("="*70)

    # Example FASTA reading
    def read_fasta(filename):
        sequence = []
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence.append(line.strip().upper())
        return ''.join(sequence)

    # Uncomment to use with actual file:
    # text = read_fasta('synthetic_datasets/SYNTH_100K/SYNTH_100K.fna')
    # pattern = "ATCGATCGATCG"
    # matches, stats = search_dna_hybrid(text, pattern, verbose=True)

    print("To use: Uncomment the code above and provide a FASTA file path")


if __name__ == "__main__":
    example_1_basic_usage()
    example_2_with_mutations()
    example_3_threshold_comparison()
    # example_4_read_fasta()  # Uncomment when you have FASTA files
