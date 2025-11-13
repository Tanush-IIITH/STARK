# KMP Algorithm: Detailed Explanation and Complexity Analysis

This document explains the Knuth–Morris–Pratt (KMP) exact string matching algorithm in detail, with special focus on the LPS (Longest Proper Prefix which is also Suffix) array construction. It also provides tight complexity analysis, edge cases, and practical notes for DNA sequence workloads.

## What KMP Solves

Given a text T of length n and a pattern P of length m, KMP finds all occurrences of P in T in O(n + m) time by avoiding re-comparison of characters that are already known to match.

- Input: text T, pattern P
- Output: all starting indices i such that T[i..i+m-1] = P
- Guarantees: linear time in the length of the text plus preprocessing of the pattern

## Notation

- n = length of the text T (n = len(T))
- m = length of the pattern P (m = len(P))

This notation is used consistently throughout the implementation and documentation.

## Core Idea

KMP precomputes, for every position i of the pattern, how far the pattern can “fallback” without losing information. This is captured by the LPS array:

- lps[i] = length of the longest proper prefix of P[0..i] which is also a suffix of P[0..i].
- “Proper” prefix means strict: not equal to the whole string itself.

During search, when a mismatch happens after j matched characters, KMP sets j = lps[j-1] instead of restarting from scratch, and continues without moving the text pointer back.

## LPS Construction (computeLPSArray)

We use the classic linear-time construction, as popularized by GeeksforGeeks (translated from C++ to Python here). The logic keeps a pointer len that tracks the candidate border length for the current position i.

Pseudocode (close to the implementation in `kmp.py`):

1. Initialize lps[0] = 0, len = 0, i = 1
2. While i < m:
   - If P[i] == P[len]:
     - len += 1
     - lps[i] = len
     - i += 1
   - Else:
     - If len != 0: set len = lps[len - 1] (fallback; do not advance i)
     - Else: lps[i] = 0; i += 1

Key properties:
- i never decreases; each step either increments i or decreases len to a previously computed value.
- Total time is O(m) because len can increase at most m times, and each decrease jumps to a strictly smaller previous border.(reduction in len is only possible at max 1 time for each increase which happens O(m) times thus reduction can also take place at max O(m) times )

### Worked Example: P = "ABABCABAB"

Index:  0 1 2 3 4 5 6 7 8
Chars:  A B A B C A B A B
lps:    0 0 1 2 0 1 2 3 4

Explanation highlights:
- At i=2 (A), it matches after len=1 → lps[2]=2? No, careful: previous matched border length is 1 ("A"); after match, len becomes 2 and we set lps[2]=2? Actually for this example, standard result yields the sequence above. The important invariant is maintained by the algorithm exactly as coded.

## Search Procedure

Maintain two indices: i for text, j for pattern.

- If T[i] == P[j]: advance both i and j
- If j == m: record a match at i - j; set j = lps[j - 1] and continue
- If mismatch and j > 0: set j = lps[j - 1]
- If mismatch and j == 0: advance i

No backtracking of i is required; hence search runs in O(n).

## Complexity Analysis

- Preprocessing (LPS): O(m)
- Search: O(n)
- Total: O(n + m)
- Space: O(m) for the LPS array

Tight bounds and notes:
- Worst-case still O(n + m); mismatches do not cause quadratic blowup because each i advances monotonically and j only falls back to previously computed borders.
- Best-case (e.g., early frequent mismatches) also linear because each character is compared at most a constant number of times.
- Constants are small: LPS computation uses only integer array writes and comparisons; search uses a couple of comparisons and jumps.

### Multiple Pattern Search

For independent patterns P1..Pk searched against the same text:
- Total time is O(n·k + sum m_i) if you re-scan the text per pattern.
- For large k, consider building a combined automaton (e.g., Aho–Corasick) to get O(n + total pattern length + output).

## DNA-Specific Notes

- DNA alphabet size is small (A, C, G, T, plus often N). KMP performance is stable regardless of alphabet size.
- Pre-normalize text/pattern to uppercase for consistent matching.
- If ambiguous bases beyond N are used, define a comparison policy (exact-only vs. IUPAC matching). The current implementation is exact-only.

## Edge Cases and Pitfalls

- Empty pattern: typically disallowed (this implementation raises ValueError in constructor).
- Pattern longer than text: search returns empty list.
- Highly repetitive patterns (e.g., "AAAA…"): KMP still runs in linear time; LPS captures the overlaps efficiently.
- Case sensitivity: the implementation normalizes to uppercase.

## Proof Sketch of Linearity

1. Preprocessing: len increases at most m times and each failed attempt sets len to lps[len-1], which is strictly smaller. Therefore, the while loop does O(m) total steps.
2. Search: pointer i over T increases monotonically from 0 to n-1. Pointer j over P never causes i to decrease; j only falls back to earlier borders via lps. Each character comparison advances i or decreases j to a previously seen border. Total comparisons are O(n).

## References

- Original papers: Knuth, D. E.; Morris, J. H.; Pratt, V. R. (1977). "Fast Pattern Matching in Strings". SIAM Journal on Computing.
- GeeksforGeeks: KMP Algorithm for Pattern Searching (computeLPSArray and search). Used as reference for LPS construction idea and overall KMP control flow. https://www.geeksforgeeks.org/kmp-algorithm-for-pattern-searching/

## Implementation Notes in This Repository

- File: `kmp.py`
  - Class `KMP` exposes `search`, `search_first`, `count_matches`, and `get_statistics`.
  - LPS is computed in `_compute_lps` and aliased by `computeLPSArray` for clarity with the GfG naming.
- File: `benchmark.py`
  - Contains microbenchmarks and simple tests mirroring the Boyer–Moore suite.
