"""
Levenshtein Fuzzy Search Implementation
Provides fuzzy pattern matching with edit distance threshold.
"""

def fuzzy_levenshtein_search(text: str, pattern: str, k: int) -> list[int]:
    """
    Finds all matches with at most 'k' Levenshtein edits using dynamic programming.
    
    Args:
        text: The text to search in
        pattern: The pattern to search for
        k: Maximum number of edits allowed (insertions, deletions, substitutions)
    
    Returns:
        List of match positions (0-indexed) where pattern matches with <= k edits
    
    Time Complexity: O(n*m) where n=len(text), m=len(pattern)
    Space Complexity: O(n) - optimized using single row
    """
    n = len(text)
    m = len(pattern)
    matches = []
    
    if m == 0 or n == 0:
        return matches
    
    # Initialize DP table with single row (space optimization)
    previous_row = [0] * (n + 1)
    
    for i in range(1, m + 1):
        current_row = [i] * (n + 1)
        for j in range(1, n + 1):
            deletion = previous_row[j] + 1
            insertion = current_row[j - 1] + 1
            substitution = previous_row[j - 1] + (0 if pattern[i - 1] == text[j - 1] else 1)
            current_row[j] = min(deletion, insertion, substitution)
        
        # On the last row, check for matches within k edits
        if i == m:
            for j in range(1, n + 1):
                if current_row[j] <= k:
                    matches.append(j - 1)
        
        previous_row = current_row
    
    return matches
