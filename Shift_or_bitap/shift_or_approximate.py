"""
Shift-Or (Bitap) Approximate String Matching Algorithm for DNA Sequences

This module implements the Shift-Or approximate matching algorithm that allows
up to k errors (substitutions, insertions, deletions) between pattern and text.

Limitation: Pattern length must be ≤ 64 base pairs (single machine word)

Notation:
- n = len(text)
- m = len(pattern) where m ≤ 64
- k = maximum number of errors allowed

Author: DNA Pattern Matching Project
Date: November 2025
"""

from typing import List, Tuple, Dict
import time

class ShiftOrApproximate:
    """
    Shift-Or approximate matching algorithm implementation (≤64 bp, k errors).

    Attributes:
        pattern (str): The pattern to search for
        pattern_length (int): Length of the pattern (m)
        k (int): Maximum number of errors allowed
        bitmasks (Dict[str, int]): Character bitmasks for bit-parallel matching
        alphabet (set): Set of characters in the pattern
        bitmask_construction_time (float): Time taken to build bitmasks
    """

    def __init__(self, pattern: str, k: int = 1):
        """
        Initialize the Shift-Or approximate matching algorithm.

        Args:
            pattern: The DNA pattern to search for (must be ≤ 64 bp)
            k: Maximum number of errors allowed (1, 2, or 3)

        Raises:
            ValueError: If pattern length > 64 or k not in [1,2,3]
        """
        if not pattern:
            raise ValueError("Pattern cannot be empty")

        if len(pattern) > 64:
            raise ValueError(f"Pattern length {len(pattern)} exceeds 64 bp limit")

        if k not in [1, 2, 3]:
            raise ValueError(f"k must be 1, 2, or 3 (got {k})")

        # Validate DNA sequence
        valid_chars = set('ACGTN')
        if not all(c in valid_chars for c in pattern.upper()):
            raise ValueError("Pattern contains invalid DNA characters")

        self.pattern = pattern.upper()
        self.pattern_length = len(pattern)
        self.k = k
        self.alphabet = set(self.pattern)
        self.bitmasks = {}
        self.bitmask_construction_time = 0.0

        # Build bitmasks
        self._build_bitmasks()

    def _build_bitmasks(self):
        """
        Build character bitmasks for the pattern.
        Same as exact matching.
        """
        start_time = time.perf_counter()

        default_mask = (1 << self.pattern_length) - 1

        for char in 'ACGTN':
            self.bitmasks[char] = default_mask

        for i, char in enumerate(self.pattern):
            self.bitmasks[char] &= ~(1 << i)

        self.bitmask_construction_time = (time.perf_counter() - start_time) * 1000

    def search(self, text: str) -> List[Tuple[int, int]]:
        """
        Search for approximate matches allowing up to k errors.

        Args:
            text: The DNA text to search in

        Returns:
            List of tuples (position, error_level) where pattern matches

        Algorithm:
        - Maintain k+1 state vectors: D0, D1, ..., Dk
        - D0 for exact matches (0 errors)
        - Di for matches with at most i errors
        - Update rules handle substitutions, insertions, deletions

        Time Complexity: O(k * n)
        Space Complexity: O(k)
        """
        if not text:
            return []

        text = text.upper()
        matches = []

        # Initialize state vectors: D0, D1, ..., Dk
        D = [(1 << self.pattern_length) - 1 for _ in range(self.k + 1)]

        match_bit = 1 << (self.pattern_length - 1)

        for j in range(len(text)):
            char = text[j]
            B_char = self.bitmasks.get(char, (1 << self.pattern_length) - 1)

            # Save old D values for computing next states
            old_D = D.copy()

            # Update D0 (exact matching)
            D[0] = ((old_D[0] << 1) | 1) & B_char

            # Update D1, D2, ..., Dk (approximate matching)
            for d in range(1, self.k + 1):
                # Substitution: old_D[d-1] << 1
                substitution = (old_D[d - 1] << 1)

                # Insertion: old_D[d]
                insertion = old_D[d]

                # Deletion: old_D[d-1]
                deletion = old_D[d - 1]

                # Match: old_D[d] << 1
                match = (old_D[d] << 1)

                # Combine all operations
                D[d] = ((substitution | insertion | deletion | match) | 1) & B_char

            # Check for matches at each error level
            for d in range(self.k + 1):
                if (D[d] & match_bit) == 0:
                    match_pos = j - self.pattern_length + 1
                    matches.append((match_pos, d))
                    break  # Report only the lowest error level

        return matches

    def search_with_metrics(self, text: str) -> Dict:
        """
        Search with detailed metrics tracking.

        Args:
            text: The DNA text to search in

        Returns:
            Dictionary with matches, time, and bit operation count
        """
        if not text:
            return {
                'matches': [],
                'search_time_ms': 0.0,
                'bit_operations': 0,
                'state_vectors': self.k + 1
            }

        text = text.upper()
        matches = []
        bit_ops = 0

        start_time = time.perf_counter()

        D = [(1 << self.pattern_length) - 1 for _ in range(self.k + 1)]
        match_bit = 1 << (self.pattern_length - 1)

        for j in range(len(text)):
            char = text[j]
            B_char = self.bitmasks.get(char, (1 << self.pattern_length) - 1)

            old_D = D.copy()

            # D0 update: 1 shift, 1 OR, 1 AND = 3 ops
            D[0] = ((old_D[0] << 1) | 1) & B_char
            bit_ops += 3

            # Each Di update: 4 shifts, 3 ORs, 1 AND = 8 ops
            for d in range(1, self.k + 1):
                substitution = (old_D[d - 1] << 1)
                insertion = old_D[d]
                deletion = old_D[d - 1]
                match = (old_D[d] << 1)
                D[d] = ((substitution | insertion | deletion | match) | 1) & B_char
                bit_ops += 8

            # Check for matches
            for d in range(self.k + 1):
                bit_ops += 1  # AND operation
                if (D[d] & match_bit) == 0:
                    match_pos = j - self.pattern_length + 1
                    matches.append((match_pos, d))
                    break

        search_time = (time.perf_counter() - start_time) * 1000

        return {
            'matches': matches,
            'search_time_ms': search_time,
            'bit_operations': bit_ops,
            'state_vectors': self.k + 1
        }

    def get_preprocessing_time(self) -> float:
        """Return bitmask construction time in milliseconds."""
        return self.bitmask_construction_time

    def get_space_usage(self) -> int:
        """
        Estimate space usage in bytes.

        Returns:
            Approximate memory usage in bytes
        """
        # Bitmasks: |alphabet| * 8 bytes
        # Pattern string: m bytes
        # State vectors: (k+1) * 8 bytes
        # Other overhead: ~100 bytes
        return len(self.bitmasks) * 8 + len(self.pattern) + (self.k + 1) * 8 + 100


def search_multiple_patterns(text: str, patterns: List[str], k: int = 1) -> Dict[str, List[Tuple[int, int]]]:
    """
    Search for multiple patterns with approximate matching.

    Args:
        text: The DNA text to search in
        patterns: List of DNA patterns to search for
        k: Maximum number of errors allowed

    Returns:
        Dictionary mapping each pattern to its list of (position, error_level) tuples
    """
    results = {}

    for pattern in patterns:
        try:
            matcher = ShiftOrApproximate(pattern, k=k)
            matches = matcher.search(text)
            results[pattern] = matches
        except ValueError as e:
            results[pattern] = []

    return results
