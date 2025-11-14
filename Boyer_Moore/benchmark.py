"""Benchmarking for Boyer-Moore on real DNA datasets.

This script:
- Discovers a genome FASTA/FNA file under STARK/DnA_dataset/ncbi_dataset/data
- Runs two analyses and writes a CSV for downstream plotting:
  1) Varying text lengths (n): fixed pattern, measure construction time/memory and search time
  2) Varying pattern lengths (m): fixed text, measure construction time/memory and search time

It is inspired by the benchmarking structure used in the Suffix Arrays-Trees module.
"""

from __future__ import annotations

import csv
import glob
import os
import random
import time
from typing import Dict, List, Tuple

import psutil

# Local imports from this folder
from boyer_moore import BoyerMoore
from utils import read_fasta_file


# Output CSV will be saved in the same directory as this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_CSV_FILE = os.path.join(SCRIPT_DIR, "bm_benchmark_results.csv")


def current_memory_bytes() -> int:
    """Return resident set size (RSS) for the current process."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss


def find_first_dataset_file() -> str | None:
    """Find a DNA FASTA/FNA file in DnA_dataset/ncbi_dataset/data (recursive).

    Returns the first matching file path or None.
    """
    script_dir = os.path.dirname(__file__)
    dataset_root = os.path.abspath(
        os.path.join(script_dir, "..", "DnA_dataset", "ncbi_dataset", "data")
    )

    patterns = ["**/*.fna", "**/*.fa", "**/*.fasta", "**/*genomic.fna"]
    candidates: List[str] = []
    for pat in patterns:
        candidates.extend(glob.glob(os.path.join(dataset_root, pat), recursive=True))

    # Prefer smaller files to keep benchmarks snappy
    candidates.sort(key=lambda p: os.path.getsize(p))
    return candidates[0] if candidates else None


def parse_fasta_contiguous(filepath: str) -> str:
    """Return the contiguous uppercase DNA string from a FASTA/FNA file."""
    seq_parts: List[str] = []
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip().upper())
    return "".join(seq_parts)


def bm_construct_and_measure(pattern: str) -> Tuple[BoyerMoore, float, int]:
    """Construct BoyerMoore object and return (obj, time_sec, memory_bytes_delta)."""
    start_mem = current_memory_bytes()
    start = time.perf_counter()
    bm = BoyerMoore(pattern)
    elapsed = time.perf_counter() - start
    end_mem = current_memory_bytes()
    return bm, elapsed, (end_mem - start_mem)


def bm_search_time(bm: BoyerMoore, text: str) -> float:
    """Measure time to find all occurrences in text."""
    start = time.perf_counter()
    _ = bm.search(text)
    return time.perf_counter() - start


def write_header(writer: csv.writer) -> None:
    writer.writerow([
        "algorithm",
        "dataset_name",
        "text_length_n",
        "pattern_length_m",
        "metric_type",
        "time_sec",
        "memory_bytes",
        "matches"
    ])


def benchmark_varying_text_lengths(full_text: str, dataset_name: str, writer: csv.writer) -> None:
    """Benchmark with a fixed pattern and varying text length n."""
    # Fixed biologically meaningful pattern (EcoRI site)
    fixed_pattern = "GAATTC"
    n_values = [
        1_000,
        2_000,
        5_000,
        10_000,
        25_000,
        50_000,
        75_000,
        100_000,
        150_000,
        200_000,
    ]

    # Always include full length if smaller than last bucket
    if len(full_text) < n_values[-1]:
        n_values = [len(full_text)]

    bm, c_time, c_mem = bm_construct_and_measure(fixed_pattern)
    # Record the construction once (pattern-only)
    writer.writerow([
        "boyer_moore",
        dataset_name,
        0,
        len(fixed_pattern),
        "construction",
        c_time,
        c_mem,
        0,
    ])

    for n in n_values:
        if n > len(full_text):
            continue
        text = full_text[:n]
        t_time = bm_search_time(bm, text)
        # Matches count (optional)
        matches = len(bm.search(text))
        writer.writerow([
            "boyer_moore",
            dataset_name,
            n,
            len(fixed_pattern),
            "search",
            t_time,
            0,
            matches,
        ])


def benchmark_varying_pattern_lengths(full_text: str, dataset_name: str, writer: csv.writer) -> None:
    """Benchmark with fixed text length and varying pattern length m."""
    random.seed(42)
    # Fix text length at up to 100k (or full length if smaller)
    n_fixed = min(100_000, len(full_text))
    text = full_text[:n_fixed]

    pattern_lengths = [5, 10, 20, 50, 100, 200]

    for m in pattern_lengths:
        if m >= n_fixed:
            continue
        # Choose a deterministic slice within text for reproducibility
        start_idx = (n_fixed // 3) % (n_fixed - m)
        pattern = text[start_idx:start_idx + m]

        bm, c_time, c_mem = bm_construct_and_measure(pattern)
        writer.writerow([
            "boyer_moore",
            dataset_name,
            n_fixed,
            m,
            "construction",
            c_time,
            c_mem,
            0,
        ])

        t_time = bm_search_time(bm, text)
        matches = len(bm.search(text))
        writer.writerow([
            "boyer_moore",
            dataset_name,
            n_fixed,
            m,
            "search",
            t_time,
            0,
            matches,
        ])


def main() -> None:
    dataset_file = find_first_dataset_file()
    if not dataset_file:
        print("No dataset found under DnA_dataset/ncbi_dataset/data. Please add a FASTA/FNA file.")
        return

    dataset_name = os.path.basename(dataset_file)
    print(f"Loading dataset: {dataset_name}")

    # Read the contiguous sequence
    full_text = parse_fasta_contiguous(dataset_file)
    print(f"Sequence length: {len(full_text):,} bp")

    with open(OUTPUT_CSV_FILE, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        write_header(writer)

        benchmark_varying_text_lengths(full_text, dataset_name, writer)
        benchmark_varying_pattern_lengths(full_text, dataset_name, writer)

    print(f"Benchmarks complete. Results saved to {OUTPUT_CSV_FILE}.")


if __name__ == "__main__":
    main()

