"""
Hybrid DNA Pattern Matcher - Architecture A (Heuristic-Driven)

Dynamically switches between:
- STATE 1: Boyer-Moore (fast with skipping)
- STATE 2: Shift-Or EXACT (bit-parallel exact matching)

Both algorithms perform EXACT matching only.
Switching based on Partial Match Density (PMD) heuristic.

Author: STARK_5 Analysis Group
Date: November 2025
"""

from boyer_moore import BoyerMoore
from shift_or_exact import shift_or_exact_search


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
    Heuristic-driven hybrid DNA pattern matcher (EXACT MATCHING ONLY).

    STATE 1 (Cruiser): Boyer-Moore for fast exact matching with skipping
    STATE 2 (Investigator): Shift-Or for bit-parallel exact matching

    Switching trigger: Partial Match Density (PMD) threshold

    When PMD is high (≥ threshold), switch to Shift-Or for efficient
    bit-parallel verification in regions with high similarity.
    """

    def __init__(self, pmd_threshold=0.75):
        """
        Initialize the hybrid matcher.

        Args:
            pmd_threshold: PMD threshold for triggering Shift-Or (default: 0.75)
        """
        self.pmd_threshold = pmd_threshold

        # Statistics tracking
        self.stats = {
            'boyer_moore_scans': 0,
            'shift_or_triggers': 0,
            'boyer_moore_matches': 0,
            'shift_or_matches': 0,
            'total_matches': 0,
            'total_positions_scanned': 0,
            'total_characters_compared': 0,
            'total_skips': 0
        }

    def search(self, text, pattern, verbose=False):
        """
        Search for pattern in text using hybrid approach.

        Args:
            text: DNA sequence to search in
            pattern: DNA pattern to search for
            verbose: Print detailed state transitions

        Returns:
            List of positions where exact matches occur
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
            print(f"HYBRID DNA MATCHER - Exact Matching Only")
            print(f"{'='*70}")
            print(f"Text length: {n:,} bp")
            print(f"Pattern length: {m} bp")
            print(f"PMD Threshold: {self.pmd_threshold}")
            print(f"Algorithm: Boyer-Moore ↔ Shift-Or (both exact)")
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
                # Exact match found by Boyer-Moore
                matches.append(i)
                self.stats['boyer_moore_matches'] += 1
                self.stats['total_matches'] += 1

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
                # High PMD but not exact → might be close match
                # STATE 2: Shift-Or (The Investigator)
                self.stats['shift_or_triggers'] += 1

                if verbose:
                    print(f"[POS {i:6d}] 🔍 PMD={pmd:.2f} ({matched_chars}/{total_chars}) "
                          f"→ TRIGGER Shift-Or")

                # Extract window for Shift-Or exact matching
                # We check a small window around current position
                window_start = max(0, i - 1)
                window_end = min(i + m + 1, n)
                window = text[window_start:window_end]

                # Run Shift-Or EXACT matching on window
                try:
                    shift_or_positions = shift_or_exact_search(window, pattern)

                    if shift_or_positions:
                        for local_pos in shift_or_positions:
                            global_pos = window_start + local_pos
                            if global_pos not in matches:  # Avoid duplicates
                                matches.append(global_pos)
                                self.stats['shift_or_matches'] += 1
                                self.stats['total_matches'] += 1

                                if verbose:
                                    print(f"[POS {global_pos:6d}] ✓ EXACT MATCH (Shift-Or)")
                except Exception as e:
                    if verbose:
                        print(f"[POS {i:6d}] ⚠ Shift-Or error: {e}")

                i += 1  # Move forward by 1 after Shift-Or trigger
            else:
                # Low PMD → use Boyer-Moore shift for fast skipping
                if j >= 0:
                    bad_char = text[i + j]
                    shift = bm.bad_character_shift(j, bad_char)
                else:
                    shift = 1

                self.stats['total_skips'] += (shift - 1)

                if verbose and shift > 1:
                    print(f"[POS {i:6d}] ⏩ SKIP {shift} positions "
                          f"(PMD={pmd:.2f} < {self.pmd_threshold})")

                i += shift

        if verbose:
            print(f"\n{'='*70}")
            print(f"Search Complete")
            print(f"{'='*70}")
            self.print_statistics()

        # Sort matches and remove duplicates
        matches = sorted(set(matches))
        return matches

    def print_statistics(self):
        """Print search statistics."""
        print(f"\n📊 Search Statistics:")
        print(f"  Boyer-Moore scans: {self.stats['boyer_moore_scans']:,}")
        print(f"  Shift-Or triggers: {self.stats['shift_or_triggers']:,}")
        print(f"  ")
        print(f"  Matches found by Boyer-Moore: {self.stats['boyer_moore_matches']}")
        print(f"  Matches found by Shift-Or: {self.stats['shift_or_matches']}")
        print(f"  Total matches: {self.stats['total_matches']}")
        print(f"  ")
        print(f"  Total positions scanned: {self.stats['total_positions_scanned']:,}")
        print(f"  Total characters compared: {self.stats['total_characters_compared']:,}")
        print(f"  Total positions skipped: {self.stats['total_skips']:,}")

        if self.stats['boyer_moore_scans'] > 0:
            trigger_rate = (self.stats['shift_or_triggers'] / 
                           self.stats['boyer_moore_scans'] * 100)
            print(f"  ")
            print(f"  Shift-Or trigger rate: {trigger_rate:.2f}%")

        # Calculate efficiency
        if self.stats['total_positions_scanned'] > 0:
            avg_comparisons = (self.stats['total_characters_compared'] / 
                             self.stats['total_positions_scanned'])
            print(f"  Avg comparisons per position: {avg_comparisons:.2f}")

            skip_rate = (self.stats['total_skips'] / 
                        self.stats['total_positions_scanned'] * 100)
            print(f"  Position skip rate: {skip_rate:.2f}%")

    def reset_statistics(self):
        """Reset all statistics counters."""
        for key in self.stats:
            self.stats[key] = 0

    def get_statistics(self):
        """Return statistics dictionary."""
        return self.stats.copy()


def search_dna_hybrid(text, pattern, pmd_threshold=0.75, verbose=False):
    """
    Convenience function for hybrid DNA pattern matching (EXACT only).

    Args:
        text: DNA sequence to search in
        pattern: DNA pattern to search for
        pmd_threshold: PMD threshold for triggering Shift-Or
        verbose: Print detailed output

    Returns:
        Tuple of (matches, statistics)
        matches: List of positions where exact matches occur
        statistics: Dictionary of search statistics
    """
    matcher = HybridDNAMatcher(pmd_threshold)
    matches = matcher.search(text, pattern, verbose)
    return matches, matcher.get_statistics()
