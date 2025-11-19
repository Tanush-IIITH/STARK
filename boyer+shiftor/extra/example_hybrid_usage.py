"""
Example Usage of Hybrid DNA Matcher

Demonstrates Architecture A (Heuristic-Driven Hybrid)
"""

from hybrid_dna_matcher import HybridDNAMatcher, search_dna_hybrid


def example_1_basic_usage():
    """Basic usage example."""
    print("="*70)
    print("EXAMPLE 1: Basic Hybrid Search")
    print("="*70)

    # DNA sequence
    text = "ACGTACGTACGTAAGTACGTACGTACGTACGT"
    pattern = "ACGTACGT"

    # Search with default parameters
    matches = search_dna_hybrid(text, pattern, verbose=True)

    print(f"\nFound {len(matches)} matches:")
    for pos, match_type, errors in matches:
        print(f"  Position {pos}: {match_type} (errors={errors})")


def example_2_with_mutations():
    """Example with mutations."""
    print("\n\n" + "="*70)
    print("EXAMPLE 2: Search with Mutations")
    print("="*70)

    # Text with mutations
    text = "ACGTACGGGACGTACGTAAGTACGTCCGTACGT"  # Mutations at positions 7, 26
    pattern = "ACGTACGT"

    # Create matcher with k=1 errors allowed
    matcher = HybridDNAMatcher(pmd_threshold=0.75, k_errors=1)
    matches = matcher.search(text, pattern, verbose=True)

    print(f"\nMatches found: {len(matches)}")


def example_3_real_genome_simulation():
    """Simulate search on a longer sequence."""
    print("\n\n" + "="*70)
    print("EXAMPLE 3: Simulated Genome Search")
    print("="*70)

    import random
    random.seed(42)

    # Generate a 10,000 bp random DNA sequence
    nucleotides = 'ACGT'
    text = ''.join(random.choice(nucleotides) for _ in range(10000))

    # Insert pattern at known positions
    pattern = "ATCGATCGATCG"
    text = text[:1000] + pattern + text[1000+len(pattern):]  # Insert at pos 1000
    text = text[:5000] + pattern + text[5000+len(pattern):]  # Insert at pos 5000

    # Search
    matcher = HybridDNAMatcher(pmd_threshold=0.80, k_errors=2)
    matches = matcher.search(text, pattern, verbose=False)

    print(f"\nText length: {len(text):,} bp")
    print(f"Pattern: {pattern}")
    print(f"Matches found: {len(matches)}")

    for pos, match_type, errors in matches:
        print(f"  Position {pos}: {match_type} (errors={errors})")

    matcher.print_statistics()


def example_4_compare_thresholds():
    """Compare different PMD thresholds."""
    print("\n\n" + "="*70)
    print("EXAMPLE 4: PMD Threshold Comparison")
    print("="*70)

    text = "ACGTACGGGACGTACGTAAGTACGTCCGTACGT" * 10
    pattern = "ACGTACGT"

    thresholds = [0.6, 0.75, 0.85, 0.95]

    for threshold in thresholds:
        print(f"\n--- PMD Threshold: {threshold} ---")
        matcher = HybridDNAMatcher(pmd_threshold=threshold, k_errors=1)
        matches = matcher.search(text, pattern, verbose=False)

        print(f"Matches: {len(matches)}")
        print(f"Exact: {matcher.stats['exact_matches']}, "
              f"Approx: {matcher.stats['approximate_matches']}")
        print(f"Shift-Or triggers: {matcher.stats['shift_or_triggers']}")


if __name__ == "__main__":
    example_1_basic_usage()
    example_2_with_mutations()
    example_3_real_genome_simulation()
    example_4_compare_thresholds()
