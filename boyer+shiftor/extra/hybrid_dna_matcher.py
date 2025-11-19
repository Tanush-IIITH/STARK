"""
Hybrid DNA Pattern Matcher - Architecture A (Heuristic-Driven)

Combines Boyer-Moore (fast exact matching) with Shift-Or (approximate matching)
using a heuristic trigger based on partial match density.

State 1 (Cruiser): Boyer-Moore for fast skipping
State 2 (Investigator): Shift-Or for approximate matching

Author: STARK_5 Analysis Group
Date: November 2025
"""

from boyer_moore import (
    boyer_moore_with_partial_match_tracking,
    calculate_partial_match_density
)
from shift_or_approximate import shift_or_approximate_search


class HybridDNAMatcher:
    """
    Heuristic-driven hybrid DNA pattern matcher.

    Dynamically switches between Boyer-Moore and Shift-Or based on
    partial match density (PMD) heuristic.
    """

    def __init__(self, pmd_threshold=0.75, k_errors=1):
        """
        Initialize the hybrid matcher.

        Args:
            pmd_threshold: PMD threshold for triggering Shift-Or (default: 0.75)
            k_errors: Maximum errors allowed in approximate matching (default: 1)
        """
        self.pmd_threshold = pmd_threshold
        self.k_errors = k_errors

        # Statistics tracking
        self.stats = {
            'boyer_moore_scans': 0,
            'shift_or_triggers': 0,
            'exact_matches': 0,
            'approximate_matches': 0,
            'total_positions_scanned': 0
        }

    def search(self, text, pattern, verbose=False):
        """
        Search for pattern in text using hybrid approach.

        Args:
            text: DNA sequence to search in
            pattern: DNA pattern to search for
            verbose: Print detailed state transitions

        Returns:
            List of tuples: (position, match_type, error_count)
            match_type: 'exact' or 'approximate'
        """
        n = len(text)
        m = len(pattern)

        if m > n or m == 0:
            return []

        matches = []
        i = 0  # Current position in text

        if verbose:
            print(f"\n{'='*70}")
            print(f"HYBRID DNA MATCHER - Starting Search")
            print(f"{'='*70}")
            print(f"Text length: {n:,} bp")
            print(f"Pattern length: {m} bp")
            print(f"PMD Threshold: {self.pmd_threshold}")
            print(f"Max errors (k): {self.k_errors}")
            print(f"{'='*70}\n")

        while i <= n - m:
            self.stats['total_positions_scanned'] += 1

            # STATE 1: Boyer-Moore (The Cruiser)
            self.stats['boyer_moore_scans'] += 1

            match_found, matched_chars, total_chars, shift = \
                boyer_moore_with_partial_match_tracking(text, pattern, i)

            if match_found:
                # Exact match found
                matches.append((i, 'exact', 0))
                self.stats['exact_matches'] += 1

                if verbose:
                    print(f"[POS {i:6d}] ✓ EXACT MATCH (Boyer-Moore)")

                i += 1  # Continue searching
                continue

            # Calculate Partial Match Density (PMD)
            pmd = calculate_partial_match_density(matched_chars, total_chars)

            # THE HEURISTIC TRIGGER
            if pmd >= self.pmd_threshold:
                # High PMD but not exact match → potential mutation site
                # STATE 2: Shift-Or (The Investigator)
                self.stats['shift_or_triggers'] += 1

                if verbose:
                    print(f"[POS {i:6d}] 🔍 PMD={pmd:.2f} → TRIGGER Shift-Or")

                # Extract local window for approximate matching
                window_start = i
                window_end = min(i + m + self.k_errors, n)
                window = text[window_start:window_end]

                # Run Shift-Or approximate matching
                approx_matches = shift_or_approximate_search(
                    window, 
                    pattern, 
                    self.k_errors
                )

                if approx_matches:
                    # Approximate match found
                    for local_pos, errors in approx_matches:
                        global_pos = window_start + local_pos
                        matches.append((global_pos, 'approximate', errors))
                        self.stats['approximate_matches'] += 1

                        if verbose:
                            print(f"[POS {global_pos:6d}] ≈ APPROX MATCH "
                                  f"(Shift-Or, k={errors})")

                i += 1  # Move forward by 1
            else:
                # Low PMD → skip using Boyer-Moore shift
                if verbose and shift > 1:
                    print(f"[POS {i:6d}] ⏩ SKIP {shift} positions "
                          f"(PMD={pmd:.2f} < {self.pmd_threshold})")

                i += shift

        if verbose:
            print(f"\n{'='*70}")
            print(f"Search Complete")
            print(f"{'='*70}")
            self.print_statistics()

        return matches

    def print_statistics(self):
        """Print search statistics."""
        print(f"\n📊 Search Statistics:")
        print(f"  Boyer-Moore scans: {self.stats['boyer_moore_scans']:,}")
        print(f"  Shift-Or triggers: {self.stats['shift_or_triggers']:,}")
        print(f"  Exact matches found: {self.stats['exact_matches']}")
        print(f"  Approximate matches found: {self.stats['approximate_matches']}")
        print(f"  Total matches: {self.stats['exact_matches'] + self.stats['approximate_matches']}")
        print(f"  Total positions scanned: {self.stats['total_positions_scanned']:,}")

        if self.stats['boyer_moore_scans'] > 0:
            trigger_rate = (self.stats['shift_or_triggers'] / 
                           self.stats['boyer_moore_scans'] * 100)
            print(f"  Shift-Or trigger rate: {trigger_rate:.2f}%")

    def reset_statistics(self):
        """Reset all statistics counters."""
        for key in self.stats:
            self.stats[key] = 0


def search_dna_hybrid(text, pattern, pmd_threshold=0.75, k_errors=1, verbose=False):
    """
    Convenience function for hybrid DNA pattern matching.

    Args:
        text: DNA sequence to search in
        pattern: DNA pattern to search for
        pmd_threshold: PMD threshold for triggering Shift-Or
        k_errors: Maximum errors allowed in approximate matching
        verbose: Print detailed output

    Returns:
        List of tuples: (position, match_type, error_count)
    """
    matcher = HybridDNAMatcher(pmd_threshold, k_errors)
    return matcher.search(text, pattern, verbose)
