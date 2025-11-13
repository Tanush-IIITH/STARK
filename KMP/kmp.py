"""
Knuth-Morris-Pratt (KMP) String Matching Algorithm for DNA Sequences

This module implements the KMP exact pattern matching algorithm from scratch,
optimized for DNA sequence analysis. It uses the prefix function (LPS array)
for linear-time searching.

Notation:
- n = len(text)
- m = len(pattern)

Reference and attribution:
- LPS construction and overall KMP control flow follow the classic
    computeLPSArray and search procedure as explained by GeeksforGeeks
    (ported from C++ to Python). See:
    https://www.geeksforgeeks.org/kmp-algorithm-for-pattern-searching/

Author: DNA Pattern Matching Project
Date: November 2025
"""

from typing import List, Tuple, Dict


class KMP:
    """
    KMP algorithm implementation for exact pattern matching in DNA sequences.

    Attributes:
        pattern (str): The pattern to search for
        pattern_length (int): Length of the pattern
        lps (List[int]): Longest Proper Prefix which is also Suffix array
    """

    def __init__(self, pattern: str):
        """
        Initialize the KMP algorithm with a pattern.

        Args:
            pattern (str): The pattern string to search for (DNA sequence)

        Raises:
            ValueError: If pattern is empty
        """
        if not pattern:
            raise ValueError("Pattern cannot be empty")

        # Normalize to uppercase for DNA strings
        # Let m = len(pattern)
        self.pattern = pattern.upper()
        self.pattern_length = len(self.pattern)  # m
        self.lps = self._compute_lps(self.pattern)

    def _compute_lps(self, pat: str) -> List[int]:
        """
        Compute the LPS (Longest Proper Prefix which is also Suffix) array.

        For each position i, lps[i] is the length of the longest proper prefix
        of pat[0..i] that is also a suffix of pat[0..i].

        Args:
            pat (str): The pattern

        Returns:
            List[int]: The LPS array
        """
    # Implementation adapted / ported from GeeksforGeeks:
        # "computeLPSArray()" used in the KMP algorithm
        # Source: https://www.geeksforgeeks.org/kmp-algorithm-for-pattern-searching/
        # The original C++ logic is preserved in behavior and translated to Python.

        m = len(pat)  # length of pattern
        lps = [0] * m
        # length of the previous longest prefix suffix
        len_pref = 0
        i = 1

        while i < m:
            if pat[i] == pat[len_pref]:
                len_pref += 1
                lps[i] = len_pref
                i += 1
            else:
                if len_pref != 0:
                    # fall back in the pattern (do not increment i)
                    len_pref = lps[len_pref - 1]
                else:
                    lps[i] = 0
                    i += 1

        return lps

    # Alias kept for clarity and traceability to the GeeksforGeeks implementation
    def computeLPSArray(self, pattern: str) -> List[int]:
        """
        computeLPSArray(pattern) -> List[int]

        Public alias for the internal LPS computation. This function name and
        algorithmic steps follow the GeeksforGeeks reference implementation
        (ported from C++ to Python).

        Reference:
        GeeksforGeeks - KMP Algorithm (computeLPSArray)
        https://www.geeksforgeeks.org/kmp-algorithm-for-pattern-searching/
        """
        return self._compute_lps(pattern)

    def search(self, text: str) -> List[int]:
        """
        Search for all occurrences of the pattern in the text.

        Args:
            text (str): The text to search in (DNA sequence)

        Returns:
            List[int]: List of starting positions where the pattern is found (0-indexed)
        
        Notes:
            Uses variables n = len(text) and m = len(pattern) consistently.
            Algorithmic structure follows the standard KMP search described by GfG.
        """
        text = text.upper()
        n = len(text)          # length of text
        m = self.pattern_length  # length of pattern
        if m == 0 or n < m:
            return []

        result = []
        i = 0  # index for text
        j = 0  # index for pattern

        while i < n:
            if self.pattern[j] == text[i]:
                i += 1
                j += 1

            if j == m:
                # Found a match ending at i-1, start index is i - j
                result.append(i - j)
                j = self.lps[j - 1]  # Continue searching for next match
            elif i < n and self.pattern[j] != text[i]:
                # mismatch after j matches
                if j != 0:
                    j = self.lps[j - 1]
                else:
                    i += 1
        return result

    def search_first(self, text: str) -> int:
        """
        Search for the first occurrence of the pattern in the text.

        Args:
            text (str): The text to search in (DNA sequence)

        Returns:
            int: Starting position of the first match, or -1 if not found
        
        Notes:
            Uses variables n = len(text) and m = len(pattern) consistently.
        """
        text = text.upper()
        n = len(text)          # length of text
        m = self.pattern_length  # length of pattern
        i = j = 0
        while i < n:
            if self.pattern[j] == text[i]:
                i += 1
                j += 1
                if j == m:
                    return i - j
            else:
                if j != 0:
                    j = self.lps[j - 1]
                else:
                    i += 1
        return -1

    def count_matches(self, text: str) -> int:
        """
        Count the total number of pattern occurrences in the text.

        Args:
            text (str): The text to search in (DNA sequence)

        Returns:
            int: Total number of matches found
        """
        return len(self.search(text))

    def get_statistics(self) -> Dict[str, any]:
        """
        Get statistics about the pattern and preprocessing.

        Returns:
            Dict[str, any]: Dictionary containing pattern statistics
        """
        return {
            'pattern': self.pattern,
            'pattern_length': self.pattern_length,
            'lps': self.lps,
        }


def search_multiple_patterns(text: str, patterns: List[str]) -> Dict[str, List[int]]:
    """
    Search for multiple patterns in the same text using the KMP algorithm.

    Args:
        text (str): The text to search in (DNA sequence)
        patterns (List[str]): List of patterns to search for

    Returns:
        Dict[str, List[int]]: Dictionary mapping each pattern to its match positions
    """
    results = {}
    for pattern in patterns:
        if pattern:
            kmp = KMP(pattern)
            results[pattern] = kmp.search(text)
    return results


def find_approximate_matches(text: str, pattern: str, max_mismatches: int = 0) -> List[Tuple[int, int]]:
    """
    Find matches allowing a specified number of mismatches.

    Note: KMP is an exact matching algorithm; for approximate matching we fall back
    to a simple sliding-window verification counting mismatches.

    Args:
        text (str): The text to search in (DNA sequence)
        pattern (str): The pattern to search for
        max_mismatches (int): Maximum number of mismatches allowed (default: 0)

    Returns:
        List[Tuple[int, int]]: List of tuples (position, number_of_mismatches)
    """
    text = text.upper()
    pattern = pattern.upper()
    n = len(text)      # length of text
    m = len(pattern)   # length of pattern
    if m == 0 or n < m:
        return []

    if max_mismatches == 0:
        kmp = KMP(pattern)
        return [(pos, 0) for pos in kmp.search(text)]

    matches: List[Tuple[int, int]] = []
    for i in range(n - m + 1):
        mismatches = 0
        for j in range(m):
            if text[i + j] != pattern[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            matches.append((i, mismatches))
    return matches
