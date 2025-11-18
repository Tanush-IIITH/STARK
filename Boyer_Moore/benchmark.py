"""Benchmark Boyer–Moore against Python's regex across real and synthetic genomes."""

from __future__ import annotations

import argparse
import csv
import glob
import os
import random
import re
import time
from dataclasses import dataclass
from typing import List, Sequence, Tuple

import psutil

from boyer_moore import BoyerMoore
from synthetic import DEFAULT_LENGTH as SYNTH_DEFAULT_LENGTH, generate_sequence
from utils import read_fasta_single_sequence


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUTPUT_CSV = os.path.join(SCRIPT_DIR, "bm_benchmark_results.csv")
FASTA_EXTENSIONS = (".fna", ".fa", ".fasta", ".ffn", ".faa", ".frn")

# Benchmark knobs mirror the suffix-array module for easy comparisons
TEXT_LENGTH_BUCKETS: Sequence[int] = (
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
)
PATTERN_LENGTHS: Sequence[int] = (5, 10, 20, 50, 100, 200)
FIXED_PATTERN = "GAATTC"  # EcoRI site
MAX_PREFIX_LENGTH = 250_000
N_FIXED_FOR_PATTERN_SWEEP = 100_000


@dataclass
class DatasetRecord:
    name: str
    group: str  # logical cohort label (e.g., "synthetic")
    text: str
    source_path: str | None = None


def current_memory_bytes() -> int:
    """Return resident set size (RSS) for the current process."""

    process = psutil.Process(os.getpid())
    return process.memory_info().rss


def discover_dataset_files(root: str, group: str) -> List[Tuple[str, str, str]]:
    """Return (path, dataset_name, group) triples for FASTA-like files."""

    if not os.path.isdir(root):
        return []

    paths: List[str] = []
    for ext in FASTA_EXTENSIONS:
        pattern = os.path.join(root, "**", f"*{ext}")
        paths.extend(glob.glob(pattern, recursive=True))

    unique_paths = sorted({os.path.abspath(p) for p in paths})
    return [(path, os.path.basename(path), group) for path in unique_paths]


def load_sequence(path: str, max_prefix: int) -> str:
    text = read_fasta_single_sequence(path).upper()
    return text[:max_prefix]


def ensure_synthetic_records(
    max_prefix: int, desired: int = 2, seeds: Sequence[int] = (137, 911)
) -> List[DatasetRecord]:
    """Build in-memory synthetic datasets when no files are present."""

    records: List[DatasetRecord] = []
    for idx in range(min(desired, len(seeds))):
        length = min(SYNTH_DEFAULT_LENGTH, max_prefix)
        sequence = generate_sequence(length, seed=seeds[idx])
        name = f"SYNTH_auto_len{length}_seed{seeds[idx]}"
        records.append(
            DatasetRecord(
                name=name,
                group="synthetic",
                text=sequence[:max_prefix],
                source_path=None,
            )
        )
    return records


def collect_datasets(dataset_root: str, max_files: int | None, max_prefix: int) -> List[DatasetRecord]:
    """Load datasets from Boyer_Moore/dataset (or a custom root)."""

    dataset_records: List[DatasetRecord] = []
    dataset_candidates = discover_dataset_files(dataset_root, "synthetic")

    if dataset_candidates:
        if max_files is not None:
            dataset_candidates = dataset_candidates[:max_files]
        for path, name, group in dataset_candidates:
            try:
                text = load_sequence(path, max_prefix)
            except FileNotFoundError:
                continue
            if not text:
                continue
            dataset_records.append(DatasetRecord(name=name, group=group, text=text, source_path=path))
    else:
        dataset_records.extend(ensure_synthetic_records(max_prefix=max_prefix))

    return dataset_records


def bm_construct_and_measure(pattern: str) -> Tuple[BoyerMoore, float, int]:
    start_mem = current_memory_bytes()
    start = time.perf_counter()
    bm = BoyerMoore(pattern)
    elapsed = time.perf_counter() - start
    end_mem = current_memory_bytes()
    return bm, elapsed, end_mem - start_mem


def bm_search_stats(bm: BoyerMoore, text: str) -> Tuple[float, int]:
    start = time.perf_counter()
    matches = bm.search(text)
    elapsed = time.perf_counter() - start
    return elapsed, len(matches)


def python_regex_search(pattern: str, text: str) -> Tuple[float, int]:
    start = time.perf_counter()
    matches = re.findall(pattern, text)
    elapsed = time.perf_counter() - start
    return elapsed, len(matches)


def write_header(writer: csv.writer) -> None:
    writer.writerow(
        [
            "algorithm",
            "dataset_name",
            "dataset_group",
            "text_length_n",
            "pattern_length_m",
            "metric_type",
            "time_sec",
            "memory_bytes",
            "matches",
        ]
    )


def benchmark_varying_text_lengths(record: DatasetRecord, writer: csv.writer) -> None:
    usable_text = record.text
    if not usable_text:
        return

    n_values = [n for n in TEXT_LENGTH_BUCKETS if n <= len(usable_text)]
    if not n_values:
        n_values = [len(usable_text)]
    elif len(usable_text) not in n_values and len(usable_text) <= max(TEXT_LENGTH_BUCKETS):
        n_values.append(len(usable_text))

    bm, c_time, c_mem = bm_construct_and_measure(FIXED_PATTERN)
    writer.writerow(
        [
            "boyer_moore",
            record.name,
            record.group,
            0,
            len(FIXED_PATTERN),
            "construction",
            c_time,
            c_mem,
            0,
        ]
    )

    for n in sorted(n_values):
        text = usable_text[:n]
        bm_time, bm_matches = bm_search_stats(bm, text)
        writer.writerow(
            [
                "boyer_moore",
                record.name,
                record.group,
                n,
                len(FIXED_PATTERN),
                "search",
                bm_time,
                0,
                bm_matches,
            ]
        )

        regex_time, regex_matches = python_regex_search(FIXED_PATTERN, text)
        writer.writerow(
            [
                "python_regex",
                record.name,
                record.group,
                n,
                len(FIXED_PATTERN),
                "search",
                regex_time,
                0,
                regex_matches,
            ]
        )


def benchmark_varying_pattern_lengths(record: DatasetRecord, writer: csv.writer) -> None:
    if not record.text:
        return

    n_fixed = min(N_FIXED_FOR_PATTERN_SWEEP, len(record.text))
    if n_fixed == 0:
        return

    text = record.text[:n_fixed]
    rng = random.Random(42)

    for m in PATTERN_LENGTHS:
        if m >= n_fixed:
            continue
        start_idx = rng.randrange(0, n_fixed - m)
        pattern = text[start_idx : start_idx + m]

        bm, c_time, c_mem = bm_construct_and_measure(pattern)
        writer.writerow(
            [
                "boyer_moore",
                record.name,
                record.group,
                n_fixed,
                m,
                "construction",
                c_time,
                c_mem,
                0,
            ]
        )

        bm_time, bm_matches = bm_search_stats(bm, text)
        writer.writerow(
            [
                "boyer_moore",
                record.name,
                record.group,
                n_fixed,
                m,
                "search",
                bm_time,
                0,
                bm_matches,
            ]
        )

        regex_time, regex_matches = python_regex_search(pattern, text)
        writer.writerow(
            [
                "python_regex",
                record.name,
                record.group,
                n_fixed,
                m,
                "search",
                regex_time,
                0,
                regex_matches,
            ]
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark Boyer-Moore vs Python regex on genome datasets stored under Boyer_Moore/dataset"
    )
    parser.add_argument(
        "--dataset-root",
        type=str,
        default=os.path.join(SCRIPT_DIR, "dataset"),
        help="Directory containing FASTA/FNA files to benchmark (default: Boyer_Moore/dataset)",
    )
    parser.add_argument(
        "--max-files",
        "--max-synthetic",
        dest="max_files",
        type=int,
        default=None,
        help="Limit the number of dataset files to benchmark (default: all)",
    )
    parser.add_argument(
        "--max-prefix",
        type=int,
        default=MAX_PREFIX_LENGTH,
        help="Maximum number of characters loaded per dataset",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=DEFAULT_OUTPUT_CSV,
        help="Path to the benchmark CSV output",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    datasets = collect_datasets(args.dataset_root, args.max_files, args.max_prefix)
    if not datasets:
        print("No datasets found in the specified folder. Add FASTA/FNA files or run synthetic.py to generate some.")
        return

    print(f"Discovered {len(datasets)} dataset(s) for benchmarking.")
    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        write_header(writer)

        for record in datasets:
            display_name = f"{record.name} ({record.group})"
            if record.source_path:
                rel_path = os.path.relpath(record.source_path, SCRIPT_DIR)
                print(f"\nDataset: {display_name} — using {rel_path}")
            else:
                print(f"\nDataset: {display_name} — synthetic in-memory sequence")

            print(f"  Available characters: {len(record.text):,}")
            benchmark_varying_text_lengths(record, writer)
            benchmark_varying_pattern_lengths(record, writer)

    print(f"\nBenchmarks complete. Results saved to {args.output}.")


if __name__ == "__main__":
    main()

