"""
Benchmarking and Testing Harness for Shift-Or Approximate Matching Algorithm

This module provides testing and benchmarking capabilities for the Shift-Or
approximate string matching algorithm on DNA sequences (≤64 bp, k=1,2,3 errors).

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
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from shift_or_approximate import ShiftOrApproximate, search_multiple_patterns
from shift_or_utils import read_fasta_file, read_fasta_single_sequence

class ShiftOrApproximateBenchmark:
    """
    Benchmarking suite for Shift-Or approximate matching (≤64 bp, k errors).
    """

    def __init__(self, k: int = 1):
        self.k = k
        self.results = []

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
        """Benchmark with increasing text sizes."""
        if slice_sizes is None:
            slice_sizes = [100000, 250000, 500000, 750000, 1000000]

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
                matcher = ShiftOrApproximate(pattern, k=self.k)
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
                'peak_memory_mb': avg_memory,
                'k': self.k
            })

        return results

    def benchmark_pattern_lengths(self, text: str,
                                  pattern_lengths: List[int] = None,
                                  trials: int = 10):
        """Benchmark impact of varying pattern lengths."""
        if pattern_lengths is None:
            pattern_lengths = [5, 10, 20, 50, 64]

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
                matcher = ShiftOrApproximate(pattern, k=self.k)
                result, elapsed, peak_mem = self.measure_time_and_memory(
                    matcher.search, text
                )
                times.append(elapsed)
                memories.append(peak_mem)
                if trial == 0:
                    matches_found = len(result)

            # Get metrics
            matcher = ShiftOrApproximate(pattern, k=self.k)
            metrics = matcher.search_with_metrics(text)

            results.append({
                'pattern_length': length,
                'n': len(text),
                'k': self.k,
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

    def benchmark_multiple_patterns(self, text: str, k_patterns: int = 50,
                                   base_length: int = 10, step: int = 5):
        """Benchmark searching for multiple patterns with k errors."""
        if len(text) > 1000000:
            text = text[:1000000]

        patterns = []
        for i in range(k_patterns):
            length = min(base_length + step * i, 64)
            start = i * 200
            if start + length <= len(text):
                pattern = text[start:start + length]
                patterns.append(pattern)

        tracemalloc.start()
        start_time = time.perf_counter()

        results = search_multiple_patterns(text, patterns, k=self.k)

        end_time = time.perf_counter()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        total_time_ms = (end_time - start_time) * 1000
        peak_memory_mb = peak / (1024 ** 2)
        total_matches = sum(len(matches) for matches in results.values())

        # Create details with error levels
        details = {}
        for pattern, matches in results.items():
            if matches:
                # Group by error level
                error_levels = {}
                for pos, error_level in matches:
                    if error_level not in error_levels:
                        error_levels[error_level] = 0
                    error_levels[error_level] += 1
                details[pattern] = {
                    'total': len(matches),
                    'by_error_level': error_levels
                }
            else:
                details[pattern] = {'total': 0, 'by_error_level': {}}

        return {
            'n': len(text),
            'k_errors': self.k,
            'k_patterns': len(patterns),
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
        """Generate scaling plots."""
        if not results:
            return

        slice_sizes = [r['slice_size'] for r in results]
        times = [r['avg_time_s'] for r in results]
        memories = [r['peak_memory_mb'] for r in results]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        ax1.plot(slice_sizes, times, 'b-o', label='Average Time')
        ax1.set_xlabel('Genome Slice Size (bp)')
        ax1.set_ylabel('Average Time (s)')
        ax1.set_title(f'Shift-Or Approximate (k={self.k}): Time vs Slice Size - {dataset_name}')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        ax2.plot(slice_sizes, memories, 'r-o', label='Peak Memory')
        ax2.set_xlabel('Genome Slice Size (bp)')
        ax2.set_ylabel('Peak Memory (MB)')
        ax2.set_title(f'Memory vs Slice Size (k={self.k})')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()

    def plot_pattern_length(self, results, output_path, dataset_name):
        """Generate pattern length plots with 64 bp marker."""
        if not results:
            return

        lengths = [r['pattern_length'] for r in results]
        mean_times = [r['time_seconds_mean'] * 1000 for r in results]
        min_times = [r['time_seconds_min'] * 1000 for r in results]
        max_times = [r['time_seconds_max'] * 1000 for r in results]
        memories = [r['peak_memory_mb'] for r in results]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        ax1.plot(lengths, mean_times, 'b-o', label='mean')
        ax1.fill_between(lengths, min_times, max_times, alpha=0.3, label='min-max')
        ax1.axvline(x=64, color='gray', linestyle='--', alpha=0.5, label='64 bp limit')
        ax1.set_xlabel('Pattern length m')
        ax1.set_ylabel('Time (ms)')
        ax1.set_title(f'Shift-Or Approximate (k={self.k}): Pattern Length Impact - {dataset_name}')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        ax2.plot(lengths, memories, 'r-o', label='Peak Memory')
        ax2.axvline(x=64, color='gray', linestyle='--', alpha=0.5, label='64 bp limit')
        ax2.set_xlabel('Pattern length m')
        ax2.set_ylabel('Peak Memory (MB)')
        ax2.set_title(f'Peak Memory vs m (k={self.k})')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()


def run_benchmark_on_dataset(fasta_file, output_dir, k_values=[1, 2, 3]):
    """
    Run complete benchmark suite for approximate matching.

    Args:
        fasta_file: Path to FASTA file
        output_dir: Directory to save results
        k_values: List of k values to test
    """
    os.makedirs(output_dir, exist_ok=True)

    sequences = read_fasta_file(fasta_file)
    if not sequences:
        return

    text = list(sequences.values())[0]
    dataset_name = os.path.basename(output_dir)

    print(f"Processing {dataset_name}...")

    for k in k_values:
        print(f"  Running benchmarks for k={k}...")
        benchmark = ShiftOrApproximateBenchmark(k=k)

        # Scaling
        scaling_results = benchmark.benchmark_scaling(text, pattern_length=20, num_runs=3)
        benchmark.save_csv(scaling_results, 
                          os.path.join(output_dir, f'scaling_results_k{k}.csv'))
        benchmark.plot_scaling(scaling_results,
                              os.path.join(output_dir, f'scaling_time_memory_k{k}.jpg'),
                              dataset_name)

        # Pattern lengths
        pattern_results = benchmark.benchmark_pattern_lengths(text, trials=10)
        benchmark.save_csv(pattern_results,
                          os.path.join(output_dir, f'pattern_length_results_k{k}.csv'))
        benchmark.plot_pattern_length(pattern_results,
                                     os.path.join(output_dir, f'pattern_length_time_memory_k{k}.jpg'),
                                     dataset_name)

        # Multiple patterns
        multi_results = benchmark.benchmark_multiple_patterns(text, k_patterns=50)
        benchmark.save_results(multi_results,
                              os.path.join(output_dir, f'multiple_patterns_summary_k{k}.json'))

    print(f"  ✓ Completed {dataset_name}")


if __name__ == "__main__":
    if len(sys.argv) > 2:
        fasta_file = sys.argv[1]
        output_dir = sys.argv[2]
        run_benchmark_on_dataset(fasta_file, output_dir)
    else:
        print("Usage: python benchmark_approximate.py <fasta_file> <output_dir>")
