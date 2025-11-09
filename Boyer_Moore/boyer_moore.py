"""
Boyer-Moore String Matching Algorithm for DNA Sequences

This module implements the Boyer-Moore exact pattern matching algorithm from scratch,
optimized for DNA sequence analysis. The implementation includes both the Bad Character
Rule and the Good Suffix Rule for efficient string matching.

Author: DNA Pattern Matching Project
Date: November 2025
"""

from typing import List, Tuple, Dict


class BoyerMoore:
    """
    Boyer-Moore algorithm implementation for exact pattern matching in DNA sequences.
    
    The algorithm uses two main heuristics:
    1. Bad Character Rule: Skip alignments based on mismatched characters
    2. Good Suffix Rule: Skip alignments based on matched suffixes
    
    Attributes:
        pattern (str): The pattern to search for
        pattern_length (int): Length of the pattern
        bad_char_table (Dict[str, List[int]]): Preprocessed bad character table
        good_suffix_table (List[int]): Preprocessed good suffix shift table
        border_position (List[int]): Helper array for good suffix preprocessing
    """
    
    def __init__(self, pattern: str):
        """
        Initialize the Boyer-Moore algorithm with a pattern.
        
        Args:
            pattern (str): The pattern string to search for (DNA sequence)
        
        Raises:
            ValueError: If pattern is empty
        """
        if not pattern:
            raise ValueError("Pattern cannot be empty")
        
        self.pattern = pattern.upper()  # Normalize to uppercase
        self.pattern_length = len(pattern)
        
        # Preprocess the pattern
        self.bad_char_table = self._preprocess_bad_character()
        self.good_suffix_table = [0] * self.pattern_length
        self.border_position = [0] * (self.pattern_length + 1)
        self._preprocess_good_suffix()
    
    def _preprocess_bad_character(self) -> Dict[str, List[int]]:
        """
        Preprocess the bad character rule table.
        
        For each character in the alphabet, store the rightmost occurrence
        of that character in the pattern at each position.
        
        Returns:
            Dict[str, List[int]]: Dictionary mapping characters to their rightmost
                                  occurrence positions at each index
        """
        # Initialize table for DNA alphabet (A, C, G, T, N)
        alphabet = ['A', 'C', 'G', 'T', 'N']
        bad_char_table = {char: [-1] * self.pattern_length for char in alphabet}
        
        # Fill the table
        for i in range(self.pattern_length):
            for char in alphabet:
                if i > 0:
                    bad_char_table[char][i] = bad_char_table[char][i - 1]
            
            # Update the current character's position
            current_char = self.pattern[i]
            if current_char in bad_char_table:
                bad_char_table[current_char][i] = i
        
        return bad_char_table
    
    def _preprocess_strong_suffix(self) -> None:
        """
        Preprocess for the strong good suffix rule.
        
        This computes the shift based on the longest suffix of pattern[i+1..m]
        that is also a suffix of the pattern and has a different preceding character.
        """
        m = self.pattern_length
        i = m
        j = m + 1
        self.border_position[i] = j
        
        while i > 0:
            # If character at position i-1 is not matching with j-1
            while j <= m and self.pattern[i - 1] != self.pattern[j - 1]:
                # Update good_suffix_table if not set yet
                if self.good_suffix_table[j - 1] == 0:
                    self.good_suffix_table[j - 1] = j - i
                j = self.border_position[j]
            
            i -= 1
            j -= 1
            self.border_position[i] = j
    
    def _preprocess_prefix_suffix(self) -> None:
        """
        Preprocess for case when pattern suffix is a prefix of the pattern.
        
        This handles the case where no proper suffix of pattern[i+1..m] matches
        a suffix of the pattern, but a prefix of the pattern matches.
        """
        m = self.pattern_length
        j = self.border_position[0]
        
        for i in range(m):
            # If good_suffix_table is not set for position i
            if self.good_suffix_table[i] == 0:
                self.good_suffix_table[i] = j
            
            # Update j when we reach the border position
            if i == j - 1:
                j = self.border_position[j]
    
    def _preprocess_good_suffix(self) -> None:
        """
        Preprocess the good suffix rule.
        
        Combines strong suffix and prefix-suffix preprocessing to create
        the complete good suffix shift table.
        """
        self._preprocess_strong_suffix()
        self._preprocess_prefix_suffix()
    
    def _get_bad_char_shift(self, text_char: str, pattern_index: int) -> int:
        """
        Calculate the shift based on the bad character rule.
        
        Args:
            text_char (str): The mismatched character from the text
            pattern_index (int): Current position in the pattern
        
        Returns:
            int: Number of positions to shift the pattern
        """
        text_char = text_char.upper()
        
        # Handle characters not in our alphabet
        if text_char not in self.bad_char_table:
            # Treat unknown characters as 'N'
            if 'N' in self.bad_char_table:
                text_char = 'N'
            else:
                return pattern_index + 1
        
        # Get the rightmost occurrence of the character before pattern_index
        rightmost_occurrence = self.bad_char_table[text_char][pattern_index]
        
        # Shift is the distance from current position to rightmost occurrence
        return pattern_index - rightmost_occurrence
    
    def search(self, text: str) -> List[int]:
        """
        Search for all occurrences of the pattern in the text.
        
        Args:
            text (str): The text to search in (DNA sequence)
        
        Returns:
            List[int]: List of starting positions where the pattern is found
                      (0-indexed)
        """
        text = text.upper()  # Normalize to uppercase
        text_length = len(text)
        matches = []
        
        # Start from the beginning of the text
        shift = 0
        
        while shift <= (text_length - self.pattern_length):
            # Start matching from right to left
            j = self.pattern_length - 1
            
            # Keep matching characters while they are equal
            while j >= 0 and self.pattern[j] == text[shift + j]:
                j -= 1
            
            # If pattern is found (j becomes -1)
            if j < 0:
                matches.append(shift)
                
                # Shift to find next occurrence
                # Use good suffix rule for shifting after a match
                shift += self.good_suffix_table[0]
            else:
                # Calculate shifts using both rules and take the maximum
                bad_char_shift = self._get_bad_char_shift(text[shift + j], j)
                good_suffix_shift = self.good_suffix_table[j]
                
                # Shift by the maximum to ensure we don't miss any matches
                shift += max(bad_char_shift, good_suffix_shift)
        
        return matches
    
    def search_first(self, text: str) -> int:
        """
        Search for the first occurrence of the pattern in the text.
        
        Args:
            text (str): The text to search in (DNA sequence)
        
        Returns:
            int: Starting position of the first match, or -1 if not found
        """
        text = text.upper()
        text_length = len(text)
        shift = 0
        
        while shift <= (text_length - self.pattern_length):
            j = self.pattern_length - 1
            
            while j >= 0 and self.pattern[j] == text[shift + j]:
                j -= 1
            
            if j < 0:
                return shift
            else:
                bad_char_shift = self._get_bad_char_shift(text[shift + j], j)
                good_suffix_shift = self.good_suffix_table[j]
                shift += max(bad_char_shift, good_suffix_shift)
        
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
            'alphabet_size': len(self.bad_char_table),
            'good_suffix_table': self.good_suffix_table,
            'bad_char_table_keys': list(self.bad_char_table.keys())
        }


def search_multiple_patterns(text: str, patterns: List[str]) -> Dict[str, List[int]]:
    """
    Search for multiple patterns in the same text using Boyer-Moore algorithm.
    
    Args:
        text (str): The text to search in (DNA sequence)
        patterns (List[str]): List of patterns to search for
    
    Returns:
        Dict[str, List[int]]: Dictionary mapping each pattern to its match positions
    """
    results = {}
    
    for pattern in patterns:
        if pattern:  # Skip empty patterns
            bm = BoyerMoore(pattern)
            matches = bm.search(text)
            results[pattern] = matches
    
    return results


def find_approximate_matches(text: str, pattern: str, max_mismatches: int = 0) -> List[Tuple[int, int]]:
    """
    Find matches allowing a specified number of mismatches.
    
    Note: This is an extension of Boyer-Moore for approximate matching.
    It uses Boyer-Moore to find candidate positions and then verifies with mismatches.
    
    Args:
        text (str): The text to search in (DNA sequence)
        pattern (str): The pattern to search for
        max_mismatches (int): Maximum number of mismatches allowed (default: 0)
    
    Returns:
        List[Tuple[int, int]]: List of tuples (position, number_of_mismatches)
    """
    if max_mismatches == 0:
        # Exact matching - use standard Boyer-Moore
        bm = BoyerMoore(pattern)
        matches = bm.search(text)
        return [(pos, 0) for pos in matches]
    
    text = text.upper()
    pattern = pattern.upper()
    text_length = len(text)
    pattern_length = len(pattern)
    matches = []
    
    # Slide window and count mismatches at each position
    for i in range(text_length - pattern_length + 1):
        mismatches = 0
        for j in range(pattern_length):
            if text[i + j] != pattern[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        
        if mismatches <= max_mismatches:
            matches.append((i, mismatches))
    
    return matches
