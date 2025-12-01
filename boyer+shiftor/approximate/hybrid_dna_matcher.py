"""
Hybrid DNA Pattern Matcher - Architecture A (Heuristic-Driven)

Integrates with existing BoyerMoore class and Shift-Or implementation.
Dynamically switches between Boyer-Moore and Shift-Or based on PMD heuristic.

Author: STARK_5 Analysis Group
Date: November 2025
"""

from boyer_moore import BoyerMoore
from shift_or_approximate import shift_or_approximate_search


def calculate_partial_match_density(text, pattern, position):
    """
    Calculate PMD by comparing pattern with text at given position.

    PMD = (Matched characters) / (Pattern length)

    Args:
        text: DNA sequence
        pattern: Pattern to match
        position: Starting position in text

    Returns:
        Tuple of (pmd_score, matched_chars, total_chars)
    """
    m = len(pattern)
    n = len(text)

    if position + m > n:
        return 0.0, 0, m

    # Compare from right to left (Boyer-Moore style)
    matched = 0
    for i in range(m - 1, -1, -1):
        if text[position + i] == pattern[i]:
            matched += 1
        else:
            break  # Stop at first mismatch

    pmd = matched / m if m > 0 else 0.0
    return pmd, matched, m


class HybridDNAMatcher:
    """
    Heuristic-driven hybrid DNA pattern matcher.

    STATE 1 (Cruiser): Boyer-Moore for fast exact matching
    STATE 2 (Investigator): Shift-Or for approximate matching

    Switching trigger: Partial Match Density (PMD) threshold
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
            'total_positions_scanned': 0,
            'total_characters_compared': 0
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

        # Initialize Boyer-Moore for this pattern
        bm = BoyerMoore(pattern)

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
            self.stats['boyer_moore_scans'] += 1

            # STATE 1: Boyer-Moore (The Cruiser)
            # Check for exact match at current position
            exact_match = True
            j = m - 1  # Start from end of pattern

            while j >= 0:
                self.stats['total_characters_compared'] += 1
                if pattern[j] != text[i + j]:
                    exact_match = False
                    break
                j -= 1

            if exact_match:
                # Exact match found
                matches.append((i, 'exact', 0))
                self.stats['exact_matches'] += 1

                if verbose:
                    print(f"[POS {i:6d}] ✓ EXACT MATCH (Boyer-Moore)")

                i += 1
                continue

            # No exact match - calculate PMD
            pmd, matched_chars, total_chars = calculate_partial_match_density(
                text, pattern, i
            )

            # THE HEURISTIC TRIGGER
            if pmd >= self.pmd_threshold:
                # High PMD but not exact → potential mutation site
                # STATE 2: Shift-Or (The Investigator)
                self.stats['shift_or_triggers'] += 1

                if verbose:
                    print(f"[POS {i:6d}] 🔍 PMD={pmd:.2f} ({matched_chars}/{total_chars}) "
                          f"→ TRIGGER Shift-Or")

                # Extract window for approximate matching
                window_start = i
                window_end = min(i + m + self.k_errors, n)
                window = text[window_start:window_end]

                # Run Shift-Or approximate matching
                try:
                    approx_matches = shift_or_approximate_search(
                        window, 
                        pattern, 
                        self.k_errors
                    )

                    if approx_matches:
                        for local_pos, errors in approx_matches:
                            global_pos = window_start + local_pos
                            matches.append((global_pos, 'approximate', errors))
                            self.stats['approximate_matches'] += 1

                            if verbose:
                                print(f"[POS {global_pos:6d}] ≈ APPROX MATCH "
                                      f"(Shift-Or, k={errors})")
                except Exception as e:
                    if verbose:
                        print(f"[POS {i:6d}] ⚠ Shift-Or error: {e}")

                i += 1  # Move forward by 1
            else:
                # Low PMD → use Boyer-Moore shift
                # Use bad character rule to determine shift
                if j >= 0:
                    bad_char = text[i + j]
                    shift = bm.bad_character_shift(j, bad_char)
                else:
                    shift = 1

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

    def search_all_overlapping(self, text, pattern, verbose=False):
        """
        Find all matches including overlapping ones.

        Args:
            text: DNA sequence to search in
            pattern: DNA pattern to search for
            verbose: Print detailed output

        Returns:
            List of tuples: (position, match_type, error_count)
        """
        # For overlapping matches, we always advance by 1
        # regardless of Boyer-Moore shift
        n = len(text)
        m = len(pattern)

        if m > n or m == 0:
            return []

        bm = BoyerMoore(pattern)
        matches = []

        for i in range(n - m + 1):
            self.stats['total_positions_scanned'] += 1

            # Check exact match
            if bm.match_at_position(text, i):
                matches.append((i, 'exact', 0))
                self.stats['exact_matches'] += 1
                continue

            # Calculate PMD
            pmd, _, _ = calculate_partial_match_density(text, pattern, i)

            # Trigger Shift-Or if PMD is high
            if pmd >= self.pmd_threshold:
                self.stats['shift_or_triggers'] += 1

                window_start = i
                window_end = min(i + m + self.k_errors, n)
                window = text[window_start:window_end]

                try:
                    approx_matches = shift_or_approximate_search(
                        window, pattern, self.k_errors
                    )

                    if approx_matches:
                        for local_pos, errors in approx_matches:
                            global_pos = window_start + local_pos
                            matches.append((global_pos, 'approximate', errors))
                            self.stats['approximate_matches'] += 1
                except:
                    pass

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
        print(f"  Total characters compared: {self.stats['total_characters_compared']:,}")

        if self.stats['boyer_moore_scans'] > 0:
            trigger_rate = (self.stats['shift_or_triggers'] / 
                           self.stats['boyer_moore_scans'] * 100)
            print(f"  Shift-Or trigger rate: {trigger_rate:.2f}%")

        # Calculate efficiency
        if self.stats['total_positions_scanned'] > 0:
            avg_comparisons = (self.stats['total_characters_compared'] / 
                             self.stats['total_positions_scanned'])
            print(f"  Avg comparisons per position: {avg_comparisons:.2f}")

    def reset_statistics(self):
        """Reset all statistics counters."""
        for key in self.stats:
            self.stats[key] = 0

    def get_statistics(self):
        """Return statistics dictionary."""
        return self.stats.copy()


def search_dna_hybrid(text, pattern, pmd_threshold=0.75, k_errors=1, 
                     overlapping=False, verbose=False):
    """
    Convenience function for hybrid DNA pattern matching.

    Args:
        text: DNA sequence to search in
        pattern: DNA pattern to search for
        pmd_threshold: PMD threshold for triggering Shift-Or
        k_errors: Maximum errors allowed in approximate matching
        overlapping: Find all overlapping matches (slower but complete)
        verbose: Print detailed output

    Returns:
        Tuple of (matches, statistics)
        matches: List of tuples (position, match_type, error_count)
        statistics: Dictionary of search statistics
    """
    matcher = HybridDNAMatcher(pmd_threshold, k_errors)

    if overlapping:
        matches = matcher.search_all_overlapping(text, pattern, verbose)
    else:
        matches = matcher.search(text, pattern, verbose)

    return matches, matcher.get_statistics()
