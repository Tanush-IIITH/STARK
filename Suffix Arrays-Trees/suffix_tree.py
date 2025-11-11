"""Naive suffix tree construction.

The implementation sticks to the textbook strategy: for every suffix beginning
at position *i* we start at the root, create edges when they are missing, and
advance one character at a time.  The resulting tree uses one node per prefix of
the input string, which makes it easy to inspect but expensive to build.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, List


@dataclass
class SuffixTreeNode:
    """Trie node used by the naive suffix tree."""

    children: Dict[str, "SuffixTreeNode"] = field(default_factory=dict)
    indices: List[int] = field(default_factory=list)

    def add_suffix(self, text: str, start_index: int) -> None:
        """Insert the suffix ``text[start_index:]`` starting from this node."""

        current = self
        for offset in range(start_index, len(text)):
            char = text[offset]
            if char not in current.children:
                current.children[char] = SuffixTreeNode()
            current = current.children[char]
            current.indices.append(start_index)


def build_naive_suffix_tree(text: str) -> SuffixTreeNode:
    """Construct the naive suffix tree for *text*.

    The tree is stored as an explicit trie where every edge represents a single
    character.  Each node records the list of suffix starting indices that reach
    it, which helps with debugging and with simple pattern matching routines.
    """

    root = SuffixTreeNode()
    for start_index in range(len(text)):
        root.add_suffix(text, start_index)

    # Storing every suffix origin at the root helps resolve empty-pattern queries.
    root.indices = list(range(len(text)))

    return root


def naive_search(root: SuffixTreeNode, pattern: str) -> List[int]:
    """Return all starting indices for *pattern* using the naive suffix tree."""

    current = root
    for char in pattern:
        if char not in current.children:
            return []
        current = current.children[char]

    return sorted(current.indices)


def _interactive_cli() -> None:
    """Drive a small REPL for building and querying the naive suffix tree."""

    text = input("Enter the text: ")
    pattern = input("Enter the pattern to search: ")

    start_time = time.perf_counter()
    tree = build_naive_suffix_tree(text)
    elapsed = time.perf_counter() - start_time

    print()
    print(f"Suffix tree built in {elapsed:.6f} seconds")

    matches = naive_search(tree, pattern)
    print()
    print(f"Pattern occurrences ({len(matches)} matches): {matches}")


if __name__ == "__main__":
    try:
        _interactive_cli()
    except KeyboardInterrupt:
        print("\nInterrupted by user.")


__all__ = [
    "SuffixTreeNode",
    "build_naive_suffix_tree",
    "naive_search",
]
