"""
Benchmark Python "re" Module vs Shift-Or/Bitap for DNA Pattern Matching

- Runs exact pattern matching using the Python re module
- Results can be compared directly with bit-parallel benchmarks

Author: Shift-Or/Bitap Analysis
Date: November 2025
"""

import os
import re
import time
import tracemalloc
import pandas as pd
from pathlib import Path

# Configuration
DATASET_ROOT = "synthetic_datasets"  # Change to "genomic_datasets" or other folder for real genomes
OUTPUT_ROOT = "re_benchmark_results"
os.makedirs(OUTPUT_ROOT, exist_ok=True)

PATTERN_LENGTHS = [5, 10, 20, 32, 50, 64, 80, 128, 256, 400]  # Same as used for Bitap
N_REPEATS = 3


def read_sequence_from_fasta(fasta_path):
    seq = []
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)


def time_python_re(sequence, pattern):
    regex = re.compile(pattern)
    tracemalloc.start()
    t0 = time.perf_counter()
    match = regex.search(sequence)
    elapsed = time.perf_counter() - t0
    peak = tracemalloc.get_traced_memory()[1] / 1024 / 1024  # MB
    tracemalloc.stop()
    return elapsed, peak, bool(match)


def run_benchmarks_on_dataset(dataset_name, fasta_file, output_dir):
    print(f"Processing {dataset_name} ...")
    sequence = read_sequence_from_fasta(fasta_file)
    n = len(sequence)
    print(f"  Sequence length: {n:,} bp")

    results = []
    for pattern_len in PATTERN_LENGTHS:
        # Extract first pattern_len basepairs as test pattern
        pattern = sequence[:pattern_len]
        if "N" in pattern:
            pattern = pattern.replace("N", "A")  # Handle ambiguous base (make valid)
        # Use plain string as regex (no wildcards)
        pattern_regex = re.escape(pattern)
        print(f"    Pattern length: {pattern_len} ...", end=" ")

        times = []
        memories = []
        matches = []
        for _ in range(N_REPEATS):
            t, m, found = time_python_re(sequence, pattern_regex)
            times.append(t)
            memories.append(m)
            matches.append(found)

        avg_time = sum(times) / N_REPEATS
        std_time = (sum((x - avg_time) ** 2 for x in times) / N_REPEATS) ** 0.5
        avg_mem = sum(memories) / N_REPEATS
        found_any = any(matches)
        print(f"mean_time={avg_time*1000:.2f} ms; mem={avg_mem:.2f} MB; found={found_any}")

        results.append({
            "dataset_name": dataset_name,
            "pattern_length": pattern_len,
            "n": n,
            "avg_time_s": avg_time,
            "std_time_s": std_time,
            "peak_memory_mb": avg_mem,
            "found": found_any
        })

    # Save result CSV
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(output_dir, "pattern_length_results_re.csv"), index=False)
    print(f"    ✓ Saved: {output_dir}/pattern_length_results_re.csv")
    


def main():
    print("="*70)
    print("PYTHON RE MODULE BENCHMARKING (EXACT MATCH)")
    print("="*70)
    print(f"Input datasets folder: {DATASET_ROOT}")
    print(f"Output: {OUTPUT_ROOT}/[dataset]/pattern_length_results_re.csv")
    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    datasets = [d for d in os.listdir(DATASET_ROOT) if os.path.isdir(os.path.join(DATASET_ROOT, d))]
    for dataset_name in sorted(datasets):
        fasta_file = None
        for file in os.listdir(os.path.join(DATASET_ROOT, dataset_name)):
            if file.endswith('.fna'):
                fasta_file = os.path.join(DATASET_ROOT, dataset_name, file)
                break
        if not fasta_file:
            print(f"  Skipping {dataset_name} - no FASTA file")
            continue
        output_dir = os.path.join(OUTPUT_ROOT, dataset_name)
        os.makedirs(output_dir, exist_ok=True)
        run_benchmarks_on_dataset(dataset_name, fasta_file, output_dir)
    print("All Python re benchmarks complete!")


if __name__ == "__main__":
    main()
