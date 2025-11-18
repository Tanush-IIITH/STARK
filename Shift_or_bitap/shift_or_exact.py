"""
Shift-Or (Bitap) Exact String Matching Algorithm for DNA Sequences

This module implements the Shift-Or (also known as Bitap) exact pattern matching 
algorithm from scratch, optimized for DNA sequence analysis. It uses bit-parallel 
operations for efficient pattern matching.

Limitation: Pattern length must be ≤ 64 base pairs (single machine word)

Notation:
- n = len(text)
- m = len(pattern) where m ≤ 64

Author: DNA Pattern Matching Project
Date: November 2025
"""

from typing import List, Tuple, Dict
import time

class ShiftOrExact:
    """
    Shift-Or exact matching algorithm implementation for DNA sequences (≤64 bp).

    Attributes:
        pattern (str): The pattern to search for
        pattern_length (int): Length of the pattern (m)
        bitmasks (Dict[str, int]): Character bitmasks for bit-parallel matching
        alphabet (set): Set of characters in the pattern
        bitmask_construction_time (float): Time taken to build bitmasks
    """

    def __init__(self, pattern: str):
        """
        Initialize the Shift-Or algorithm with a pattern.

        Args:
            pattern: The DNA pattern to search for (must be ≤ 64 bp)

        Raises:
            ValueError: If pattern length > 64 or contains invalid characters
        """
        if not pattern:
            raise ValueError("Pattern cannot be empty")

        if len(pattern) > 64:
            raise ValueError(f"Pattern length {len(pattern)} exceeds 64 bp limit for standard Shift-Or")

        # Validate DNA sequence
        valid_chars = set('ACGTN')
        if not all(c in valid_chars for c in pattern.upper()):
            raise ValueError("Pattern contains invalid DNA characters")

        self.pattern = pattern.upper()
        self.pattern_length = len(pattern)
        self.alphabet = set(self.pattern)
        self.bitmasks = {}
        self.bitmask_construction_time = 0.0

        # Build bitmasks
        self._build_bitmasks()

    def _build_bitmasks(self):
        """
        Build character bitmasks for the pattern.

        For each character c in the alphabet:
        - B[c] has bit i set to 0 if pattern[i] == c, else 1
        - This allows efficient bit-parallel matching

        Time Complexity: O(m * |Σ|) where |Σ| is alphabet size
        Space Complexity: O(|Σ|)
        """
        start_time = time.perf_counter()

        # Initialize all bitmasks to all 1s
        default_mask = (1 << self.pattern_length) - 1

        # For DNA, we typically have A, C, G, T, N
        for char in 'ACGTN':
            self.bitmasks[char] = default_mask

        # Set bits to 0 where character matches pattern position
        for i, char in enumerate(self.pattern):
            # Clear bit i for character char
            self.bitmasks[char] &= ~(1 << i)

        self.bitmask_construction_time = (time.perf_counter() - start_time) * 1000  # ms

    def search(self, text: str) -> List[int]:
        """
        Search for all occurrences of the pattern in the text.

        Args:
            text: The DNA text to search in

        Returns:
            List of starting positions where pattern is found

        Algorithm:
        - Maintain state vector D (all 1s initially)
        - For each character, update: D = ((D << 1) | 1) & B[char]
        - Match when bit (m-1) is 0

        Time Complexity: O(n)
        Space Complexity: O(1)
        """
        if not text:
            return []

        text = text.upper()
        matches = []

        # Initialize state: all bits set to 1
        D = (1 << self.pattern_length) - 1

        # Bit position to check for match
        match_bit = 1 << (self.pattern_length - 1)

        for j in range(len(text)):
            char = text[j]

            # Get bitmask for current character (default to all 1s if not in alphabet)
            B_char = self.bitmasks.get(char, (1 << self.pattern_length) - 1)

            # Update state
            D = ((D << 1) | 1) & B_char

            # Check for match (bit m-1 is 0)
            if (D & match_bit) == 0:
                match_pos = j - self.pattern_length + 1
                matches.append(match_pos)

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
                'state_vectors': 1
            }

        text = text.upper()
        matches = []
        bit_ops = 0

        start_time = time.perf_counter()

        D = (1 << self.pattern_length) - 1
        match_bit = 1 << (self.pattern_length - 1)

        for j in range(len(text)):
            char = text[j]
            B_char = self.bitmasks.get(char, (1 << self.pattern_length) - 1)

            # Count operations: 1 shift, 1 OR, 1 AND
            D = ((D << 1) | 1) & B_char
            bit_ops += 3

            # Count operation: 1 AND
            if (D & match_bit) == 0:
                bit_ops += 1
                match_pos = j - self.pattern_length + 1
                matches.append(match_pos)

        search_time = (time.perf_counter() - start_time) * 1000  # ms

        return {
            'matches': matches,
            'search_time_ms': search_time,
            'bit_operations': bit_ops,
            'state_vectors': 1
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
        # Bitmasks: |alphabet| * 8 bytes per int
        # Pattern string: m bytes
        # Other overhead: ~100 bytes
        return len(self.bitmasks) * 8 + len(self.pattern) + 100


def search_multiple_patterns(text: str, patterns: List[str]) -> Dict[str, List[int]]:
    """
    Search for multiple patterns in the text.

    Args:
        text: The DNA text to search in
        patterns: List of DNA patterns to search for

    Returns:
        Dictionary mapping each pattern to its list of match positions
    """
    results = {}

    for pattern in patterns:
        try:
            matcher = ShiftOrExact(pattern)
            matches = matcher.search(text)
            results[pattern] = matches
        except ValueError as e:
            # Skip invalid patterns
            results[pattern] = []

    return results
