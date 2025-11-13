"""
Example usage for KMP algorithm (KMP/example_usage.py)

Shows simple examples of using the KMP class for DNA sequence pattern matching.
"""

import os
import sys

# Make sure the package directory is importable
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from kmp import KMP, search_multiple_patterns
from utils import generate_random_dna, get_reverse_complement, calculate_gc_content


def basic_example():
    text = "ACGTACGTATGCATGCACGTATGCTAGC"
    pattern = "ATGC"

    print("Text:", text)
    print("Pattern:", pattern)

    k = KMP(pattern)
    matches = k.search(text)
    print("Matches at:", matches)


def performance_example():
    text = generate_random_dna(100000, seed=42)
    pattern = text[1000:1020]
    k = KMP(pattern)

    import time
    start = time.time()
    matches = k.search(text)
    end = time.time()

    print(f"Searched text length {len(text)} in {(end-start)*1000:.3f} ms; matches: {len(matches)}")


def main():
    print("\nKMP Example Usage\n" + "="*40)
    basic_example()
    performance_example()


if __name__ == '__main__':
    main()
