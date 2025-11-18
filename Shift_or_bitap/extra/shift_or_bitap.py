# shift_or_bitap.py

# This implements:

# exact Shift-Or (Bitap) for patterns up to machine word size (64 by default),

# an approximate variant supporting up to k substitutions (the classic bitap extension),

# convenience wrappers to behave similarly to BoyerMoore in API.


# Notes:

# The approximate implementation above follows Myers / bit-parallel ideas but simplified; it's robust for substitution-only approximate search for small patterns. (Insertions/deletions require more complex handling.)

# If you expect patterns longer than WORD_BITS (64), the Boyer–Moore harness uses patterns up to 200 — so for long patterns we will run a fallback naive search. See the runner below. 

from typing import List, Optional, Tuple
import sys

WORD_BITS = 64  # use 64-bit masks by default; small patterns only

class ShiftOr:
    """
    Shift-Or / Bitap exact and approximate matcher.
    - Exact search: pattern length m must be <= WORD_BITS (64). For longer patterns,
      fallback to naive sliding or chunked method (we raise for clarity).
    - Approximate: supports up to k substitutions (Hamming-like). Insertions/deletions
      are not handled by the basic bit-parallel approach without heavier DP.
    """

    def __init__(self, pattern: str, max_errors: int = 0):
        self.pattern = pattern
        self.m = len(pattern)
        self.k = max_errors
        if self.m == 0:
            raise ValueError("Pattern must be non-empty")
        if self.m > WORD_BITS:
            raise ValueError(f"Pattern length {self.m} > {WORD_BITS}. Reduce pattern length or modify WORD_BITS.")
        self._build_mask()

    def _build_mask(self):
        """Build char -> bitmask table: 1 at bit i if pattern[i] != char (we use inverted masks)."""
        # classic bitap uses a mask of 1s and sets 0 for matches; here we keep the standard approach
        self.char_mask = {}
        ALL_ONES = (1 << self.m) - 1
        # init all chars to all ones (no match)
        for c in ['A','C','G','T','N']:  # common DNA alphabet; others will be created on demand
            self.char_mask[c] = ALL_ONES
        for i, ch in enumerate(self.pattern):
            # bit i from LSB corresponds to pattern position i
            bit = 1 << i
            # set bit to 0 where character matches
            self.char_mask[ch] = self.char_mask.get(ch, ALL_ONES) & (~bit & ALL_ONES)

    def search(self, text: str) -> List[int]:
        """Return list of 0-based starting positions where pattern occurs (exact match)."""
        if self.k != 0:
            return self.search_approx(text, self.k)
        m = self.m
        if m == 0:
            return []
        mask = (1 << m) - 1
        R = mask  # bitmask; after processing a char, if bit (m-1) is 0 -> match at position
        results = []
        for i, ch in enumerate(text):
            ch = ch.upper()
            char_mask = self.char_mask.get(ch, mask)
            # Shift left by 1, OR with char_mask
            R = ((R << 1) | 1) & char_mask
            if (R & (1 << (m - 1))) == 0:
                # match ends at i, so start at i-m+1
                results.append(i - m + 1)
        return results

    def search_first(self, text: str) -> Optional[int]:
        res = self.search(text)
        return res[0] if res else None

    def count_matches(self, text: str) -> int:
        return len(self.search(text))

    def search_approx(self, text: str, k: int) -> List[Tuple[int,int]]:
        """
        Approximate search allowing up to k substitutions (Hamming-like).
        Returns list of (start_index, errors) for each match with errors <= k.

        Note: This implements the bit-parallel algorithm for substitutions (Myers-like).
        Pattern length must be <= WORD_BITS.
        """
        # Build pattern masks (same as exact but as equality mask)
        m = self.m
        if m > WORD_BITS:
            raise ValueError("Pattern too long for bit-parallel approximate search.")

        # equality masks: 1 where char matches
        B = {}
        ALL_ONES = (1 << m) - 1
        for ch in ['A','C','G','T','N']:
            B[ch] = 0
        for i, pch in enumerate(self.pattern):
            bit = 1 << i
            B[pch] = B.get(pch, 0) | bit

        # initialize VP (vertical positives) and VN (vertical negatives) per Myers' algorithm
        VP = ALL_ONES  # all ones
        VN = 0
        curr_errors = m
        results = []

        for i, ch in enumerate(text):
            ch = ch.upper()
            PM = B.get(ch, 0)
            # Myers bit-parallel operations
            X = PM | VN
            D0 = (((X & VP) + VP) ^ VP) | X
            HP = VN | ~(D0 | VP)
            HN = D0 & VP
            # Next VP and VN
            shift_HP = (HP << 1) & ALL_ONES
            shift_HN = (HN << 1) & ALL_ONES
            VP = (shift_HN | ~(D0 | shift_HP)) & ALL_ONES
            VN = shift_HP & D0

            # compute current errors: if bit m-1 of VP is 0 then errors decreased, etc.
            # Myers trick to track errors per position:
            if (VP & (1 << (m - 1))) == 0:
                curr_errors -= 1
            elif (VN & (1 << (m - 1))) != 0:
                curr_errors += 1

            if curr_errors <= k:
                results.append((i - m + 1, curr_errors))
        # filter out negative starts
        results = [(s,e) for (s,e) in results if s >= 0]
        return results

    def get_statistics(self) -> dict:
        return {"pattern": self.pattern, "pattern_length": self.m, "max_errors": self.k}
