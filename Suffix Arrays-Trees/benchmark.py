"""Benchmark harness for suffix array and suffix tree implementations.

The script follows the project roadmap:
1. Discover FASTA datasets.
2. Build search structures (suffix arrays / suffix trees) while tracking
   runtime and memory usage.
3. Run pattern queries against the built structures alongside Python's
   regex engine for comparison.
4. Persist all measurements in ``benchmark_results.csv`` for downstream
   analysis.
"""

from __future__ import annotations

import csv
import glob
import os
import re
import time
from typing import Callable, Dict, Iterable, List, Tuple

import psutil

from suffix_array import (
    brute_force_suffix_array,
    locate_pattern,
    manber_myers_suffix_array,
)
from suffix_tree import (
    build_naive_suffix_tree,
    naive_search,
    build_ukkonen_suffix_tree,
    ukkonen_search,
)


def parse_fasta_file(filepath: str) -> str:
    """Return a contiguous uppercase string extracted from a FASTA file.

    FASTA records mix header (``>``-prefixed) lines with raw sequence lines.
    Only the sequence lines contribute to the text benchmarked by the
    algorithms, so we filter headers and concatenate the remainder.
    """

    sequence_parts: List[str] = []
    with open(filepath, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            sequence_parts.append(line.strip().upper())
    return "".join(sequence_parts)


def find_datasets() -> List[str]:
    """Discover candidate FASTA datasets relative to the project root."""

    script_dir = os.path.dirname(__file__)
    dataset_root = os.path.join(script_dir, "dataset")
    pattern = os.path.join(dataset_root, "**", "*.fna")
    return glob.glob(pattern, recursive=True)


def current_memory_bytes() -> int:
    """Return resident set size (RSS) for the current process."""

    process = psutil.Process(os.getpid())
    return process.memory_info().rss


def benchmark_construction(build_fn: Callable[[str], object], text: str) -> Tuple[object, float, int]:
    """Time and measure memory for a construction routine."""

    start_mem = current_memory_bytes()
    start_time = time.perf_counter()
    structure = build_fn(text)
    end_time = time.perf_counter()
    end_mem = current_memory_bytes()

    return structure, end_time - start_time, end_mem - start_mem


def benchmark_query(query_fn: Callable[..., object], *args: object) -> float:
    """Return execution time for a query routine."""

    start_time = time.perf_counter()
    query_fn(*args)
    end_time = time.perf_counter()
    return end_time - start_time


def write_csv_row(writer: csv.writer, row: Iterable[object]) -> None:
    """Write a single row to the CSV output."""

    writer.writerow(row)


if __name__ == "__main__":
    OUTPUT_CSV_FILE = "benchmark_results.csv"
    PATTERNS_TO_SEARCH = [
    "TCTGTGT",
    "TAAAATTTTATTGACTTA",
    "CACCACCATCACCATTACCAC",
    "TTTTTTTTTT",
    "AGCTTTTCATTCTGACTGCAACGGGCAATA",
    "GTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTT",
    "AACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGA",
    "GGCGGCAATATCGAAACTGTTGCCATCGACGGCGATTTCGATGCCTGTCAGGCGCTGGTGAAGCAGGCGTTTGATGATGAAGAACTGAAAGTGGCGCTAGGGTTAAACTCGGCTAACTCGATTAACATCAGCCGTTTGCTGGCGCAGA",
    "GATTACAZZZ",
    ]
    # Run sa_naive up to n=10000 (will run on 1k, 5k, 10k)
    SA_NAIVE_THRESHOLD: int = 100_000

    # Run st_naive only up to n=2000 (will run on 1k only)
    ST_NAIVE_THRESHOLD: int = 2_001
    N_SIZES_FOR_LARGE_FILES = [
        1_000,
        2_000,
        3_000,
        4_000,
        5_000,
        10_000,
        25_000,
        50_000,
        75_000,
        100_000,
        150_000,
        200_000,
    ]

    datasets = find_datasets()
    if not datasets:
        print("No datasets found. Ensure .fna files are located under dataset/ relative to this script.")

    with open(OUTPUT_CSV_FILE, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        write_csv_row(
            writer,
            [
                "algorithm",
                "dataset_name",
                "text_length_n",
                "pattern_length_m",
                "metric_type",
                "time_sec",
                "memory_bytes",
            ],
        )

        for dataset_path in datasets:
            dataset_name = os.path.basename(dataset_path)
            print(f"\nLoading dataset: {dataset_name}...")
            full_text = parse_fasta_file(dataset_path)
            full_length = len(full_text)

            if full_length < N_SIZES_FOR_LARGE_FILES[-1]:
                n_values_to_test = [full_length]
            else:
                n_values_to_test = N_SIZES_FOR_LARGE_FILES

            for n in n_values_to_test:
                if n > full_length:
                    continue

                text = full_text[:n]
                print(f"  Benchmarking {dataset_name} (n={n})...")

                built_structures: Dict[str, object] = {}

                if n <= SA_NAIVE_THRESHOLD:
                    structure, elapsed, memory = benchmark_construction(brute_force_suffix_array, text)
                    built_structures["sa_naive"] = structure
                    write_csv_row(
                        writer,
                        [
                            "sa_naive",
                            dataset_name,
                            n,
                            0,
                            "construction",
                            elapsed,
                            memory,
                        ],
                    )
                else:
                    print(f"    Skipping sa_naive (n > {SA_NAIVE_THRESHOLD})")

                if n <= ST_NAIVE_THRESHOLD:
                    structure, elapsed, memory = benchmark_construction(build_naive_suffix_tree, text)
                    built_structures["st_naive"] = structure
                    write_csv_row(
                        writer,
                        [
                            "st_naive",
                            dataset_name,
                            n,
                            0,
                            "construction",
                            elapsed,
                            memory,
                        ],
                    )
                else:
                    print(f"    Skipping st_naive (n > {ST_NAIVE_THRESHOLD})")

                structure, elapsed, memory = benchmark_construction(build_ukkonen_suffix_tree, text)
                built_structures["st_ukkonen"] = structure
                write_csv_row(
                    writer,
                    [
                        "st_ukkonen",
                        dataset_name,
                        n,
                        0,
                        "construction",
                        elapsed,
                        memory,
                    ],
                )

                structure, elapsed, memory = benchmark_construction(manber_myers_suffix_array, text)
                built_structures["sa_manber"] = structure
                write_csv_row(
                    writer,
                    [
                        "sa_manber_myers",
                        dataset_name,
                        n,
                        0,
                        "construction",
                        elapsed,
                        memory,
                    ],
                )

                for pattern in PATTERNS_TO_SEARCH:
                    m = len(pattern)

                    if "sa_manber" in built_structures:
                        elapsed = benchmark_query(locate_pattern, text, pattern, built_structures["sa_manber"])
                        write_csv_row(
                            writer,
                            [
                                "sa_manber_query",
                                dataset_name,
                                n,
                                m,
                                "query",
                                elapsed,
                                0,
                            ],
                        )

                    if "sa_naive" in built_structures:
                        elapsed = benchmark_query(locate_pattern, text, pattern, built_structures["sa_naive"])
                        write_csv_row(
                            writer,
                            [
                                "sa_naive_query",
                                dataset_name,
                                n,
                                m,
                                "query",
                                elapsed,
                                0,
                            ],
                        )

                    if "st_naive" in built_structures:
                        elapsed = benchmark_query(naive_search, built_structures["st_naive"], pattern)
                        write_csv_row(
                            writer,
                            [
                                "st_naive_query",
                                dataset_name,
                                n,
                                m,
                                "query",
                                elapsed,
                                0,
                            ],
                        )

                    if "st_ukkonen" in built_structures:
                        elapsed = benchmark_query(ukkonen_search, built_structures["st_ukkonen"], pattern)
                        write_csv_row(
                            writer,
                            [
                                "st_ukkonen_query",
                                dataset_name,
                                n,
                                m,
                                "query",
                                elapsed,
                                0,
                            ],
                        )

                    elapsed = benchmark_query(re.findall, pattern, text)
                    write_csv_row(
                        writer,
                        [
                            "python_regex",
                            dataset_name,
                            n,
                            m,
                            "query",
                            elapsed,
                            0,
                        ],
                    )

    print(f"Benchmarks complete. Results saved to {OUTPUT_CSV_FILE}.")
