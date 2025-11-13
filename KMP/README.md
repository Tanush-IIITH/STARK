# KMP Algorithm for DNA Pattern Matching

A from-scratch implementation of the Knuth–Morris–Pratt (KMP) exact string matching algorithm optimized for DNA sequence analysis, plus a simple benchmarking harness.

## Overview

The **KMP algorithm** achieves linear-time exact pattern matching by preprocessing the pattern into an LPS (Longest Proper Prefix which is also Suffix) array. This allows efficient backtracking-free scanning of the text.

## Notation

- n = length of the text (len(text))
- m = length of the pattern (len(pattern))

This notation is used consistently in code comments, docstrings, and this README.

## Features

- ✅ Pure Python KMP implementation (no dependencies)
- ✅ Optimized for DNA sequences (A, C, G, T, N)
- ✅ Multiple pattern search helper
- ✅ Benchmark suite mirroring the Boyer–Moore directory
- ✅ FASTA file utilities bundled for easy use

## Project Structure

```
KMP/
├── kmp.py                 # Core KMP algorithm implementation
├── utils.py               # DNA sequence utilities (FASTA I/O, random DNA, etc.)
├── benchmark.py           # Tests and benchmarks (random + optional real dataset)
└── README.md              # This file
```

## Requirements

- Python 3.7+
- No external dependencies (standard library only)

## Quick Start

Run the benchmark suite (tests + microbenchmarks):

```bash
python benchmark.py
```

This will:
- Run correctness tests and edge cases
- Benchmark by pattern length and text length
- Benchmark multiple-pattern search
- Optionally benchmark on the real dataset if found locally

A `benchmark_results.json` file will be written in the `KMP/` folder.

## Algorithm Explanation and Complexity

For a deep dive into how KMP works, including the LPS array construction, proofs of linearity, and edge cases, see:

- ALGORITHM_EXPLANATION.md (this folder)

The LPS routine is based on the classic computeLPSArray flow popularized by GeeksforGeeks (see References below), translated to Python and documented in code (`KMP.computeLPSArray`).

## Using the KMP Class

```python
from kmp import KMP

text = "ACGTACGTATGCATGCACGT"
pattern = "ATGC"

kmp = KMP(pattern)
positions = kmp.search(text)      # [8, 12]
first = kmp.search_first(text)    # 8
count = kmp.count_matches(text)   # 2
stats = kmp.get_statistics()      # {'pattern': 'ATGC', 'pattern_length': 4, 'lps': [...]}
```

Multiple pattern search helper:

```python
def search_multiple_patterns(text, patterns):
    # returns {pattern: [positions, ...], ...}
```

## Working with FASTA files

`utils.py` provides simple helpers to load sequences:

- `read_fasta_file(path)` → dict of {header: sequence}
- `read_fasta_single_sequence(path)` → concatenated string
- `read_fasta_sequences_only(path)` → list of sequences
- `get_all_fasta_files(dir)` → list of .fna/.fa/.fasta files

## Real Dataset Benchmark (optional)

If the repository also contains the dataset at:
```
STARK/DnA_dataset/ncbi_dataset/data
```
then `benchmark.py` will auto-detect it and run the real-data benchmark. If not found, that step is skipped.

You can also pass a custom path by editing `dataset_path_candidates` in `benchmark.py`.

## Notes

- KMP performs especially well for long texts and is robust for all alphabets; in DNA (small alphabet) it’s reliable and predictable.
- For approximate matching (with mismatches), `kmp.py` includes a simple fallback that scans with a mismatch counter; this is not KMP proper but handy for small allowances.

## License

Educational use for algorithms coursework.
 
## References

- KMP algorithm and LPS computation adapted from GeeksforGeeks (ported from C++):
    https://www.geeksforgeeks.org/kmp-algorithm-for-pattern-searching/
