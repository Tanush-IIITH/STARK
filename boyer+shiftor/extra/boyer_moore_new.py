"""
Boyer-Moore Algorithm for DNA Pattern Matching

Implements the Boyer-Moore algorithm with Bad Character Rule
for fast exact pattern matching.

Author: STARK_5 Analysis Group
Date: November 2025
"""

def build_bad_character_table(pattern):
    """
    Build the Bad Character Rule table for Boyer-Moore.

    Args:
        pattern: The pattern string to search for

    Returns:
        Dictionary mapping each character to its rightmost occurrence position
    """
    bad_char = {}
    m = len(pattern)

    # For DNA, we only care about A, C, G, T
    for char in 'ACGT':
        bad_char[char] = -1

    # Record rightmost occurrence of each character in pattern
    for i in range(m):
        bad_char[pattern[i]] = i

    return bad_char


def boyer_moore_search(text, pattern):
    """
    Perform Boyer-Moore exact pattern matching.

    Args:
        text: The text to search in
        pattern: The pattern to search for

    Returns:
        List of positions where exact matches occur
    """
    n = len(text)
    m = len(pattern)

    if m > n:
        return []

    # Preprocess
    bad_char = build_bad_character_table(pattern)

    matches = []
    i = 0  # Position in text

    while i <= n - m:
        j = m - 1  # Start from end of pattern

        # Compare pattern from right to left
        while j >= 0 and pattern[j] == text[i + j]:
            j -= 1

        if j < 0:
            # Exact match found
            matches.append(i)
            i += 1
        else:
            # Mismatch - use bad character rule
            bad_char_shift = j - bad_char.get(text[i + j], -1)
            i += max(1, bad_char_shift)

    return matches


def boyer_moore_with_partial_match_tracking(text, pattern, start_pos):
    """
    Enhanced Boyer-Moore that tracks partial match quality.

    Args:
        text: The text to search in
        pattern: The pattern to search for
        start_pos: Starting position in text

    Returns:
        Tuple of (match_found, matched_chars, total_chars, shift_amount)
    """
    m = len(pattern)
    n = len(text)

    if start_pos + m > n:
        return False, 0, m, 1

    # Build bad character table
    bad_char = build_bad_character_table(pattern)

    # Scan from right to left
    matched_chars = 0
    j = m - 1

    while j >= 0:
        if pattern[j] == text[start_pos + j]:
            matched_chars += 1
            j -= 1
        else:
            # Mismatch found
            break

    if j < 0:
        # Exact match
        return True, m, m, 1

    # Calculate shift using bad character rule
    bad_char_shift = j - bad_char.get(text[start_pos + j], -1)
    shift = max(1, bad_char_shift)

    return False, matched_chars, m, shift


def calculate_partial_match_density(matched_chars, total_chars):
    """
    Calculate the Partial Match Density (PMD).

    PMD = Matched Characters / Total Pattern Length

    Args:
        matched_chars: Number of characters that matched
        total_chars: Total pattern length

    Returns:
        PMD score between 0 and 1
    """
    if total_chars == 0:
        return 0.0
    return matched_chars / total_chars
