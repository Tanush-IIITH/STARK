"""
Benchmarking and Testing Harness for Shift-Or Exact Matching Algorithm

This module provides testing and benchmarking capabilities for the Shift-Or
exact string matching algorithm on DNA sequences (≤64 bp patterns).

Author: DNA Pattern Matching Project
Date: November 2025
"""

import time
import os
import sys
import tracemalloc
from typing import List, Dict, Tuple
import json
from datetime import datetime
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from shift_or_exact import ShiftOrExact, search_multiple_patterns
from shift_or_utils import (
    read_fasta_file, read_fasta_single_sequence,
    get_all_fasta_files, generate_random_dna
)

class ShiftOrExactBenchmark:
    """
    Benchmarking suite for the Shift-Or exact matching algorithm (≤64 bp).
    """

    def __init__(self, dataset_path: str = None):
        self.results = []
        self.dataset_path = dataset_path

    def measure_time_and_memory(self, func, *args, **kwargs):
        """Measure execution time and peak memory usage."""
        tracemalloc.start()
        start_time = time.perf_counter()

        result = func(*args, **kwargs)

        end_time = time.perf_counter()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        elapsed_time = end_time - start_time
        peak_memory_mb = peak / (1024 ** 2)

        return result, elapsed_time, peak_memory_mb

    def benchmark_scaling(self, text: str, pattern_length: int = 20, 
                         slice_sizes: List[int] = None, num_runs: int = 3):
        """
        Benchmark algorithm performance with increasing text sizes.

        Args:
            text: Full DNA sequence
            pattern_length: Length of pattern to extract
            slice_sizes: List of text slice sizes to test
            num_runs: Number of repetitions for averaging

        Returns:
            List of results for each slice size
        """
        if slice_sizes is None:
            slice_sizes = [100000, 250000, 500000, 750000, 1000000]

        # Extract pattern from fixed position
        if len(text) < 1000 + pattern_length:
            return []

        pattern = text[1000:1000 + pattern_length]
        results = []

        for size in slice_sizes:
            if size > len(text):
                continue

            text_slice = text[:size]
            times = []
            memories = []

            for _ in range(num_runs):
                matcher = ShiftOrExact(pattern)
                _, elapsed, peak_mem = self.measure_time_and_memory(
                    matcher.search, text_slice
                )
                times.append(elapsed)
                memories.append(peak_mem)

            avg_time = sum(times) / len(times)
            avg_memory = sum(memories) / len(memories)

            results.append({
                'slice_size': size,
                'avg_time_s': avg_time,
                'peak_memory_mb': avg_memory
            })

        return results

    def benchmark_pattern_lengths(self, text: str, 
                                  pattern_lengths: List[int] = None,
                                  trials: int = 10):
        """
        Benchmark impact of varying pattern lengths.

        Args:
            text: Full DNA sequence (fixed at 1M bp)
            pattern_lengths: List of pattern lengths to test (≤64 bp)
            trials: Number of independent trials for statistics

        Returns:
            List of results for each pattern length
        """
        if pattern_lengths is None:
            pattern_lengths = [5, 10, 20, 50, 64]

        # Use fixed text size
        if len(text) > 1000000:
            text = text[:1000000]

        results = []

        for length in pattern_lengths:
            if len(text) < 1000 + length:
                continue

            pattern = text[1000:1000 + length]
            times = []
            memories = []
            matches_found = 0

            for trial in range(trials):
                matcher = ShiftOrExact(pattern)
                result, elapsed, peak_mem = self.measure_time_and_memory(
                    matcher.search, text
                )
                times.append(elapsed)
                memories.append(peak_mem)
                if trial == 0:
                    matches_found = len(result)

            # Get metrics from matcher
            matcher = ShiftOrExact(pattern)
            metrics = matcher.search_with_metrics(text)

            results.append({
                'pattern_length': length,
                'n': len(text),
                'time_seconds_mean': sum(times) / len(times),
                'time_seconds_min': min(times),
                'time_seconds_max': max(times),
                'peak_memory_mb': sum(memories) / len(memories),
                'matches_found': matches_found,
                'bitmask_time_ms': matcher.get_preprocessing_time(),
                'bit_operations': metrics['bit_operations'],
                'state_vectors': metrics['state_vectors'],
                'implementation': 'single-word'
            })

        return results

    def benchmark_multiple_patterns(self, text: str, k: int = 50,
                                   base_length: int = 10, step: int = 5):
        """
        Benchmark searching for multiple patterns.

        Args:
            text: Full DNA sequence
            k: Number of patterns to test
            base_length: Starting pattern length
            step: Length increment for each pattern

        Returns:
            Summary dictionary with total time and matches
        """
        if len(text) > 1000000:
            text = text[:1000000]

        patterns = []
        for i in range(k):
            length = min(base_length + step * i, 64)  # Cap at 64 bp
            start = i * 200
            if start + length <= len(text):
                pattern = text[start:start + length]
                patterns.append(pattern)

        tracemalloc.start()
        start_time = time.perf_counter()

        results = search_multiple_patterns(text, patterns)

        end_time = time.perf_counter()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        total_time_ms = (end_time - start_time) * 1000
        peak_memory_mb = peak / (1024 ** 2)
        total_matches = sum(len(matches) for matches in results.values())

        # Create details dictionary
        details = {pattern: len(matches) for pattern, matches in results.items()}

        return {
            'n': len(text),
            'k': len(patterns),
            'total_time_ms': total_time_ms,
            'peak_memory_mb': peak_memory_mb,
            'total_matches': total_matches,
            'details': details
        }

    def save_results(self, results, filepath):
        """Save results to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(results, f, indent=2)

    def save_csv(self, results, filepath):
        """Save results to CSV file."""
        if not results:
            return

        import csv
        keys = results[0].keys()

        with open(filepath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(results)

    def plot_scaling(self, results, output_path, dataset_name):
        """Generate scaling plots (time and memory vs slice size)."""
        if not results:
            return

        slice_sizes = [r['slice_size'] for r in results]
        times = [r['avg_time_s'] for r in results]
        memories = [r['peak_memory_mb'] for r in results]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        # Time plot
        ax1.plot(slice_sizes, times, 'b-o', label='Average Time')
        ax1.set_xlabel('Genome Slice Size (bp)')
        ax1.set_ylabel('Average Time (s)')
        ax1.set_title(f'Shift-Or: Time vs Slice Size (n) - {dataset_name}')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Memory plot
        ax2.plot(slice_sizes, memories, 'r-o', label='Peak Memory')
        ax2.set_xlabel('Genome Slice Size (bp)')
        ax2.set_ylabel('Peak Memory (MB)')
        ax2.set_title('Memory vs Slice Size (n)')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()

    def plot_pattern_length(self, results, output_path, dataset_name):
        """Generate pattern length plots with 64 bp boundary marker."""
        if not results:
            return

        lengths = [r['pattern_length'] for r in results]
        mean_times = [r['time_seconds_mean'] * 1000 for r in results]  # Convert to ms
        min_times = [r['time_seconds_min'] * 1000 for r in results]
        max_times = [r['time_seconds_max'] * 1000 for r in results]
        memories = [r['peak_memory_mb'] for r in results]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        # Time plot with min-max range
        ax1.plot(lengths, mean_times, 'b-o', label='mean')
        ax1.fill_between(lengths, min_times, max_times, alpha=0.3, label='min-max')
        # Add 64 bp boundary marker
        ax1.axvline(x=64, color='gray', linestyle='--', alpha=0.5, label='64 bp limit')
        ax1.set_xlabel('Pattern length m')
        ax1.set_ylabel('Time (ms)')
        ax1.set_title(f'Shift-Or: Impact of Pattern Length (fixed n) - {dataset_name}')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Memory plot
        ax2.plot(lengths, memories, 'r-o', label='Peak Memory')
        ax2.axvline(x=64, color='gray', linestyle='--', alpha=0.5, label='64 bp limit')
        ax2.set_xlabel('Pattern length m')
        ax2.set_ylabel('Peak Memory (MB)')
        ax2.set_title('Peak Memory vs m')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()


def run_benchmark_on_dataset(fasta_file, output_dir):
    """
    Run complete benchmark suite on a single dataset.

    Args:
        fasta_file: Path to FASTA file
        output_dir: Directory to save results
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Read sequence
    sequences = read_fasta_file(fasta_file)
    if not sequences:
        return

    # Get first sequence
    text = list(sequences.values())[0]
    dataset_name = os.path.basename(output_dir)

    print(f"Processing {dataset_name}...")

    # Initialize benchmark
    benchmark = ShiftOrExactBenchmark()

    # 1. Scaling benchmark
    print("  Running scaling benchmark...")
    scaling_results = benchmark.benchmark_scaling(text, pattern_length=20, num_runs=3)
    benchmark.save_csv(scaling_results, os.path.join(output_dir, 'scaling_results.csv'))
    benchmark.plot_scaling(scaling_results, 
                          os.path.join(output_dir, 'scaling_time_memory.jpg'),
                          dataset_name)

    # 2. Pattern length benchmark
    print("  Running pattern length benchmark...")
    pattern_results = benchmark.benchmark_pattern_lengths(text, trials=10)
    benchmark.save_csv(pattern_results, 
                      os.path.join(output_dir, 'pattern_length_results.csv'))
    benchmark.plot_pattern_length(pattern_results,
                                 os.path.join(output_dir, 'pattern_length_time_memory.jpg'),
                                 dataset_name)

    # 3. Multiple patterns benchmark
    print("  Running multiple patterns benchmark...")
    multi_results = benchmark.benchmark_multiple_patterns(text, k=50)
    benchmark.save_results(multi_results,
                          os.path.join(output_dir, 'multiple_patterns_summary.json'))

    print(f"  ✓ Completed {dataset_name}")


if __name__ == "__main__":
    # Example usage
    if len(sys.argv) > 2:
        fasta_file = sys.argv[1]
        output_dir = sys.argv[2]
        run_benchmark_on_dataset(fasta_file, output_dir)
    else:
        print("Usage: python benchmark_exact.py <fasta_file> <output_dir>")
