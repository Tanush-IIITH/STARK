"""Suffix tree construction.

Includes:
1. Naive suffix tree construction - textbook strategy with O(n^2) complexity.
2. Ukkonen's algorithm - linear time construction with optimal O(n) complexity.
"""

from __future__ import annotations
import time
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional

# =============================================================================
# Naive Suffix Tree (Unaffected)
# =============================================================================

@dataclass
class SuffixTreeNode:
    """Trie node used by the naive suffix tree."""
    children: Dict[str, "SuffixTreeNode"] = field(default_factory=dict)
    indices: List[int] = field(default_factory=list)

    def add_suffix(self, text: str, start_index: int) -> None:
        current = self
        for offset in range(start_index, len(text)):
            char = text[offset]
            if char not in current.children:
                current.children[char] = SuffixTreeNode()
            current = current.children[char]
            current.indices.append(start_index)

def build_naive_suffix_tree(text: str) -> SuffixTreeNode:
    root = SuffixTreeNode()
    for start_index in range(len(text)):
        root.add_suffix(text, start_index)
    root.indices = list(range(len(text)))
    return root

def naive_search(root: SuffixTreeNode, pattern: str) -> List[int]:
    current = root
    for char in pattern:
        if char not in current.children:
            return []
        current = current.children[char]
    return sorted(current.indices)

# =============================================================================
# Ukkonen's Suffix Tree (Corrected Implementation)
# =============================================================================

class UkkonenNode:
    """Node class for a correct Ukkonen's suffix tree implementation."""
    def __init__(self, start: int, end_ptr: List[int], suffix_index: int = -1):
        self.start: int = start
        self.end: List[int] = end_ptr 
        self.children: Dict[str, UkkonenNode] = {}
        self.suffix_link: Optional[UkkonenNode] = None
        self.suffix_index: int = suffix_index

    def edge_len(self) -> int:
        return self.end[0] - self.start + 1

class UkkonenSuffixTree:
    """Correct O(n) Ukkonen's suffix tree construction."""

    def __init__(self, text: str):
        self.text: str = text
        self.n: int = len(text)
        self.global_end: List[int] = [-1]
        self.root: UkkonenNode = UkkonenNode(-1, [-1])
        self.root.suffix_link = self.root
        self.active_node: UkkonenNode = self.root
        self.active_edge: Optional[str] = None
        self.active_length: int = 0
        self.remainder: int = 0
        
        for i in range(self.n):
            self.tree_extend(i)

    def walk_down(self, node: UkkonenNode) -> bool:
        edge_len = node.edge_len()
        if self.active_length >= edge_len:
            self.active_length -= edge_len
            self.active_edge = self.text[node.start + edge_len]
            self.active_node = node
            return True
        return False

    def tree_extend(self, text_idx: int) -> None:
        char = self.text[text_idx]
        self.global_end[0] = text_idx
        self.remainder += 1
        last_new_node: Optional[UkkonenNode] = None

        while self.remainder > 0:
            if self.active_length == 0:
                self.active_edge = char

            if self.active_edge not in self.active_node.children:
                # Rule 2: Create a new leaf
                new_leaf = UkkonenNode(text_idx, self.global_end, suffix_index = text_idx - (self.remainder - 1))
                self.active_node.children[self.active_edge] = new_leaf
                if last_new_node:
                    last_new_node.suffix_link = self.active_node
                    last_new_node = None
            else:
                # Rule 3: Match/Split
                next_node = self.active_node.children[self.active_edge]
                if self.walk_down(next_node):
                    continue 
                if self.text[next_node.start + self.active_length] == char:
                    # Rule 3, Case 1
                    self.active_length += 1
                    if last_new_node:
                        last_new_node.suffix_link = self.active_node
                        last_new_node = None
                    break 
                
                # Rule 3, Case 2
                split_end = [next_node.start + self.active_length - 1]
                split_node = UkkonenNode(next_node.start, split_end, suffix_index = -1)
                self.active_node.children[self.active_edge] = split_node
                new_leaf = UkkonenNode(text_idx, self.global_end, suffix_index = text_idx - (self.remainder - 1))
                split_node.children[char] = new_leaf
                next_node.start += self.active_length
                split_node.children[self.text[next_node.start]] = next_node
                if last_new_node:
                    last_new_node.suffix_link = split_node
                last_new_node = split_node

            self.remainder -= 1
            if self.active_node == self.root and self.active_length > 0:
                self.active_length -= 1
                self.active_edge = self.text[text_idx - self.remainder + 1]
            else:
                self.active_node = self.active_node.suffix_link if self.active_node.suffix_link else self.root

    def _search_path(self, pattern: str) -> Optional[UkkonenNode]:
        """Traverses the tree for the pattern. Returns end node or None."""
        node = self.root
        pat_idx = 0
        while pat_idx < len(pattern):
            char = pattern[pat_idx]
            if char not in node.children:
                return None
            child = node.children[char]
            edge_len = child.edge_len()
            for i in range(edge_len):
                if pat_idx + i >= len(pattern):
                    return child
                if self.text[child.start + i] != pattern[pat_idx + i]:
                    return None
            pat_idx += edge_len
            node = child
        return node

    def _collect_leaves(self, node: UkkonenNode, results: List[int]) -> None:
        """Recursively collect suffix_index from all leaves under this node."""
        if node.suffix_index != -1:  # This is a leaf
            results.append(node.suffix_index)
            return
        for child in node.children.values():
            self._collect_leaves(child, results)

    # --- FUNCTION LOGIC CORRECTED ---
    def search_pattern(self, pattern: str) -> List[int]:
        """Return all starting indices for *pattern*."""
        
        # FIX 1: Empty pattern ("")
        # The naive tree (given "text$") returns all indices from 0 to n.
        # This implementation must do the same. self.n includes the '$'.
        if not pattern:
            return list(range(self.n)) # Changed from (self.n - 1)
        
        # 1. Find the node where the pattern path ends
        node = self._search_path(pattern)
        
        if node is None:
            return []  # Pattern not found
        
        # 2. Collect all leaf indices under that node
        results: List[int] = []
        self._collect_leaves(node, results)
        
        # FIX 2: Remove filtering
        # The naive tree (on "text$") includes all matching indices.
        # This implementation must also include all of them.
        # The previous filtering was incorrect for this test.
        return sorted(results)


def build_ukkonen_suffix_tree(text: str) -> UkkonenSuffixTree:
    """
    Construct suffix tree using Ukkonen's linear time algorithm.
    Appends a unique terminal '$' to the text.
    """
    if not text.endswith("$"):
        text += "$"
    
    tree = UkkonenSuffixTree(text)
    return tree


def ukkonen_search(tree: UkkonenSuffixTree, pattern: str) -> List[int]:
    """Return all starting indices for *pattern* using Ukkonen's suffix tree."""
    return tree.search_pattern(pattern)


# =============================================================================
# Main CLI (Unaffected)
# =============================================================================

def _interactive_cli() -> None:
    """Drive a small REPL for building and querying suffix trees."""

    # A dictionary mapping user's choice to the build and search functions
    algorithms = {
        "naive": (build_naive_suffix_tree, naive_search),
        "ukkonen": (build_ukkonen_suffix_tree, ukkonen_search)
    }

    print("Available algorithms: naive, ukkonen")
    while True:
        choice = input("Select algorithm: ").strip().lower()
        if choice in algorithms:
            break
        print("Unrecognised choice. Please enter 'naive' or 'ukkonen'.")

    text = input("Enter the text: ")
    pattern = input("Enter the pattern to search: ")

    build_func, search_func = algorithms[choice]

    # --- Build Step ---
    print(f"\nBuilding suffix tree using '{choice}'...")
    start_time = time.perf_counter()
    # Note: build_ukkonen adds '$', build_naive does not.
    # For a fair comparison, we add '$' for naive.
    build_text = text if choice == "ukkonen" else (text + "$")
    
    # build_ukkonen_suffix_tree is called if choice is "ukkonen"
    # build_naive_suffix_tree is called if choice is "naive"
    tree = build_func(build_text) 
    build_elapsed = time.perf_counter() - start_time
    print(f"Suffix tree built in {build_elapsed:.6f} seconds")

    # --- Search Step ---
    start_time = time.perf_counter()
    # ukkonen_search is called if choice is "ukkonen"
    # naive_search is called if choice is "naive"
    matches = search_func(tree, pattern) 
    search_elapsed = time.perf_counter() - start_time
    print(f"Search completed in {search_elapsed:.6f} seconds")

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
    "UkkonenState"
    "UkkonenNode",
    "UkkonenSuffixTree",
    "build_ukkonen_suffix_tree",
    "ukkonen_search",
]