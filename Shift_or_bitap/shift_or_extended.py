"""
Shift-Or (Bitap) Extended Algorithm for Long Patterns (>64 bp)

This module implements an extended multi-word version of Shift-Or for patterns
longer than 64 base pairs. It splits the pattern across multiple 64-bit words.

Limitation: Pattern length must be ≤ 800 base pairs

Notation:
- n = len(text)
- m = len(pattern) where 64 < m ≤ 800
- w = number of 64-bit words needed = ceil(m / 64)

Author: DNA Pattern Matching Project
Date: November 2025
"""

from typing import List, Tuple, Dict
import time
import math

class ShiftOrExtended:
    """
    Extended Shift-Or algorithm for patterns > 64 bp (exact matching only).

    Attributes:
        pattern (str): The pattern to search for
        pattern_length (int): Length of the pattern (m)
        num_words (int): Number of 64-bit words needed
        bitmasks (Dict[str, List[int]]): Character bitmasks for each word
        alphabet (set): Set of characters in the pattern
        bitmask_construction_time (float): Time taken to build bitmasks
    """

    def __init__(self, pattern: str):
        """
        Initialize the extended Shift-Or algorithm.

        Args:
            pattern: The DNA pattern to search for (64 < length ≤ 800 bp)

        Raises:
            ValueError: If pattern length ≤ 64 or > 800
        """
        if not pattern:
            raise ValueError("Pattern cannot be empty")

        if len(pattern) <= 64:
            raise ValueError(f"Pattern length {len(pattern)} ≤ 64 bp. Use ShiftOrExact instead.")

        if len(pattern) > 800:
            raise ValueError(f"Pattern length {len(pattern)} exceeds 800 bp limit")

        # Validate DNA sequence
        valid_chars = set('ACGTN')
        if not all(c in valid_chars for c in pattern.upper()):
            raise ValueError("Pattern contains invalid DNA characters")

        self.pattern = pattern.upper()
        self.pattern_length = len(pattern)
        self.num_words = math.ceil(self.pattern_length / 64)
        self.alphabet = set(self.pattern)
        self.bitmasks = {}
        self.bitmask_construction_time = 0.0

        # Build multi-word bitmasks
        self._build_bitmasks()

    def _build_bitmasks(self):
        """
        Build character bitmasks for multiple words.

        Each character gets a list of bitmasks, one per word.
        Word 0 covers pattern[0:64], word 1 covers pattern[64:128], etc.
        """
        start_time = time.perf_counter()

        # Initialize bitmasks for each character
        for char in 'ACGTN':
            self.bitmasks[char] = [(1 << 64) - 1 for _ in range(self.num_words)]

        # Set bits for each character position
        for i, char in enumerate(self.pattern):
            word_idx = i // 64
            bit_pos = i % 64

            # Clear bit at position bit_pos in word word_idx
            self.bitmasks[char][word_idx] &= ~(1 << bit_pos)

        # Handle last word (may not be full 64 bits)
        bits_in_last_word = self.pattern_length % 64
        if bits_in_last_word > 0:
            # Mask off unused bits in last word
            mask = (1 << bits_in_last_word) - 1
            for char in self.bitmasks:
                self.bitmasks[char][-1] &= mask

        self.bitmask_construction_time = (time.perf_counter() - start_time) * 1000

    def search(self, text: str) -> List[int]:
        """
        Search for pattern in text using multi-word Shift-Or.

        Args:
            text: The DNA text to search in

        Returns:
            List of starting positions where pattern is found

        Algorithm:
        - Maintain w state words: D[0], D[1], ..., D[w-1]
        - For each character, update all words with carry propagation
        - Match when final bit of last word is 0

        Time Complexity: O(w * n) where w = ceil(m/64)
        Space Complexity: O(w)
        """
        if not text:
            return []

        text = text.upper()
        matches = []

        # Initialize state words
        D = [(1 << 64) - 1 for _ in range(self.num_words)]

        # Initialize last word properly
        bits_in_last_word = self.pattern_length % 64
        if bits_in_last_word > 0:
            D[-1] = (1 << bits_in_last_word) - 1

        # Bit position to check for match in last word
        if bits_in_last_word > 0:
            match_bit = 1 << (bits_in_last_word - 1)
        else:
            match_bit = 1 << 63

        for j in range(len(text)):
            char = text[j]

            # Get bitmasks for current character
            B_char = self.bitmasks.get(char, None)
            if B_char is None:
                B_char = [(1 << 64) - 1 for _ in range(self.num_words)]
                if bits_in_last_word > 0:
                    B_char[-1] = (1 << bits_in_last_word) - 1

            # Update state words with carry propagation
            carry = 0
            for w in range(self.num_words):
                # Shift left and add carry from previous word
                shifted = (D[w] << 1) | carry

                # Compute carry for next word (bit 64)
                carry = 1 if (D[w] & (1 << 63)) == 0 else 0

                # Apply bitmask
                D[w] = (shifted | 1) & B_char[w]

            # Check for match in last word
            if (D[-1] & match_bit) == 0:
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
                'state_vectors': self.num_words
            }

        text = text.upper()
        matches = []
        bit_ops = 0

        start_time = time.perf_counter()

        D = [(1 << 64) - 1 for _ in range(self.num_words)]
        bits_in_last_word = self.pattern_length % 64
        if bits_in_last_word > 0:
            D[-1] = (1 << bits_in_last_word) - 1
            match_bit = 1 << (bits_in_last_word - 1)
        else:
            match_bit = 1 << 63

        for j in range(len(text)):
            char = text[j]
            B_char = self.bitmasks.get(char, None)
            if B_char is None:
                B_char = [(1 << 64) - 1 for _ in range(self.num_words)]
                if bits_in_last_word > 0:
                    B_char[-1] = (1 << bits_in_last_word) - 1

            carry = 0
            for w in range(self.num_words):
                shifted = (D[w] << 1) | carry
                carry = 1 if (D[w] & (1 << 63)) == 0 else 0
                D[w] = (shifted | 1) & B_char[w]
                # Operations: 2 shifts, 3 ORs, 2 ANDs = 7 ops per word
                bit_ops += 7

            bit_ops += 1  # AND for match check
            if (D[-1] & match_bit) == 0:
                match_pos = j - self.pattern_length + 1
                matches.append(match_pos)

        search_time = (time.perf_counter() - start_time) * 1000

        return {
            'matches': matches,
            'search_time_ms': search_time,
            'bit_operations': bit_ops,
            'state_vectors': self.num_words
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
        # Bitmasks: |alphabet| * num_words * 8 bytes
        # Pattern string: m bytes
        # State vectors: num_words * 8 bytes
        # Other overhead: ~100 bytes
        return len(self.bitmasks) * self.num_words * 8 + len(self.pattern) + self.num_words * 8 + 100


def search_multiple_patterns(text: str, patterns: List[str]) -> Dict[str, List[int]]:
    """
    Search for multiple extended patterns in the text.

    Args:
        text: The DNA text to search in
        patterns: List of DNA patterns to search for (>64 bp each)

    Returns:
        Dictionary mapping each pattern to its list of match positions
    """
    results = {}

    for pattern in patterns:
        try:
            matcher = ShiftOrExtended(pattern)
            matches = matcher.search(text)
            results[pattern] = matches
        except ValueError as e:
            results[pattern] = []

    return results
