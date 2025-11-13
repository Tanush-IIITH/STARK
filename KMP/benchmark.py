"""
Benchmarking and Testing Harness for KMP Algorithm

This module provides testing and benchmarking capabilities for the KMP
string matching algorithm on DNA sequences, mirroring the Boyer-Moore suite.

Author: DNA Pattern Matching Project
Date: November 2025
"""

import time
import os
import sys
from typing import List, Dict, Tuple
import json
from datetime import datetime

# Add current directory to path to import our modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from kmp import KMP, search_multiple_patterns, find_approximate_matches
from utils import (
    read_fasta_file, read_fasta_single_sequence, read_fasta_sequences_only,
    get_all_fasta_files, generate_random_dna, validate_dna_sequence, 
    get_reverse_complement, calculate_gc_content, read_fasta_generator,
    count_nucleotides
)


class KMPBenchmark:
    """
    Benchmarking suite for the KMP algorithm.

    Attributes:
        results (List[Dict]): List of benchmark results
        dataset_path (str): Path to the DNA dataset directory
    """

    def __init__(self, dataset_path: str = None):
        self.results: List[Dict] = []
        self.dataset_path = dataset_path

    def _measure_time(self, func, *args, **kwargs) -> Tuple[any, float]:
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        return result, end_time - start_time

    def test_basic_matching(self) -> Dict:
        print("\n=== Testing Basic Pattern Matching (KMP) ===")
        test_cases = [
            {'name': 'Simple match', 'text': 'ACGTACGTACGT', 'pattern': 'ACG', 'expected': [0, 4, 8]},
            {'name': 'No match', 'text': 'ACGTACGTACGT', 'pattern': 'TTT', 'expected': []},
            {'name': 'Multiple matches', 'text': 'ACGTACGTACGT', 'pattern': 'GTACG', 'expected': [2, 6]},
            {'name': 'Pattern equals text', 'text': 'ACGT', 'pattern': 'ACGT', 'expected': [0]},
            {'name': 'Overlapping matches', 'text': 'AAAAAAA', 'pattern': 'AAA', 'expected': [0, 1, 2, 3, 4]},
        ]

        passed = failed = 0
        details = []
        for test in test_cases:
            kmp = KMP(test['pattern'])
            result = kmp.search(test['text'])
            status = 'PASS' if result == test['expected'] else 'FAIL'
            passed += status == 'PASS'
            failed += status == 'FAIL'
            details.append({'test': test['name'], 'status': status, 'expected': test['expected'], 'got': result})
            print(f"  {test['name']}: {status}")
        print(f"\nBasic Tests: {passed} passed, {failed} failed")
        return {'test_type': 'basic_matching', 'passed': passed, 'failed': failed, 'details': details}

    def test_edge_cases(self) -> Dict:
        print("\n=== Testing Edge Cases (KMP) ===")
        passed = failed = 0

        # Single char pattern
        try:
            kmp = KMP('A')
            res = kmp.search('ACGTACGT')
            expected = [0, 4]
            if res == expected:
                print("  Single character pattern: PASS"); passed += 1
            else:
                print("  Single character pattern: FAIL"); failed += 1
        except Exception as e:
            print(f"  Single character pattern: FAIL ({e})"); failed += 1

        # Pattern longer than text
        try:
            kmp = KMP('ACGTACGTACGT')
            res = kmp.search('ACG')
            expected = []
            if res == expected:
                print("  Pattern longer than text: PASS"); passed += 1
            else:
                print("  Pattern longer than text: FAIL"); failed += 1
        except Exception as e:
            print(f"  Pattern longer than text: FAIL ({e})"); failed += 1

        # Case insensitivity
        try:
            kmp = KMP('acgt')
            res = kmp.search('ACGTACGT')
            expected = [0, 4]
            if res == expected:
                print("  Case insensitivity: PASS"); passed += 1
            else:
                print("  Case insensitivity: FAIL"); failed += 1
        except Exception as e:
            print(f"  Case insensitivity: FAIL ({e})"); failed += 1

        print(f"\nEdge Case Tests: {passed} passed, {failed} failed")
        return {'test_type': 'edge_cases', 'passed': passed, 'failed': failed}

    def benchmark_pattern_length(self, text_length: int = 100000) -> Dict:
        print("\n=== Benchmarking Pattern Length Impact (KMP) ===")
        text = generate_random_dna(text_length, seed=42)
        pattern_lengths = [5, 10, 20, 50, 100, 200]
        results = []
        for length in pattern_lengths:
            pattern = text[1000:1000+length]
            kmp = KMP(pattern)
            matches, t = self._measure_time(kmp.search, text)
            results.append({'pattern_length': length, 'text_length': text_length, 'matches_found': len(matches), 'time_seconds': t, 'time_ms': t*1000})
            print(f"  Pattern length {length:3d}: {t*1000:8.3f} ms ({len(matches)} matches)")
        return {'benchmark_type': 'pattern_length', 'text_length': text_length, 'results': results}

    def benchmark_text_length(self, pattern_length: int = 20) -> Dict:
        print("\n=== Benchmarking Text Length Impact (KMP) ===")
        text_lengths = [10000, 50000, 100000, 500000]
        pattern = generate_random_dna(pattern_length, seed=42)
        results = []
        for length in text_lengths:
            text = generate_random_dna(length, seed=100)
            kmp = KMP(pattern)
            matches, t = self._measure_time(kmp.search, text)
            cps = length / t if t > 0 else float('inf')
            results.append({'text_length': length, 'pattern_length': pattern_length, 'matches_found': len(matches), 'time_seconds': t, 'time_ms': t*1000, 'chars_per_second': cps})
            msg = f"  Text length {length:7d}: {t*1000:8.3f} ms " + (f"({length/t/1e6:.2f} M chars/sec)" if t > 0 else "(<0.001 ms)")
            print(msg)
        return {'benchmark_type': 'text_length', 'pattern_length': pattern_length, 'results': results}

    def benchmark_multiple_patterns(self, text_length: int = 100000, num_patterns: int = 10) -> Dict:
        print("\n=== Benchmarking Multiple Pattern Search (KMP) ===")
        text = generate_random_dna(text_length, seed=42)
        patterns: List[str] = []
        for i in range(num_patterns):
            length = 10 + i * 5
            patterns.append(text[i*100:i*100+length])
        matches_dict, t = self._measure_time(search_multiple_patterns, text, patterns)
        total_matches = sum(len(m) for m in matches_dict.values())
        print(f"  Searched {num_patterns} patterns in {t*1000:.3f} ms; total matches {total_matches}")
        pattern_results = [
            {'pattern': (p[:20] + '...') if len(p) > 20 else p, 'pattern_length': len(p), 'matches': len(matches_dict[p])}
            for p in patterns
        ]
        return {'benchmark_type': 'multiple_patterns', 'text_length': text_length, 'num_patterns': num_patterns, 'total_matches': total_matches, 'total_time_ms': t*1000, 'avg_time_per_pattern_ms': (t*1000/num_patterns if num_patterns else 0), 'pattern_results': pattern_results}

    def benchmark_real_dataset(self, max_files: int = 3) -> Dict:
        print("\n=== Benchmarking on Real DNA Dataset (KMP) ===")
        if not self.dataset_path or not os.path.exists(self.dataset_path):
            print("  Dataset path not found. Skipping real dataset benchmark.")
            return {'benchmark_type': 'real_dataset', 'status': 'skipped'}
        fasta_files = get_all_fasta_files(self.dataset_path, recursive=True)
        if not fasta_files:
            print("  No FASTA files found in dataset. Skipping.")
            return {'benchmark_type': 'real_dataset', 'status': 'no_files'}
        fasta_files = fasta_files[:max_files]
        print(f"  Processing {len(fasta_files)} files...")

        patterns = [
            'ATGCATGC', 'GCTAGCTA', 'TATAAA', 'CAAT', 'GAATTC', 'GGATCC'
        ]

        results = []
        for filepath in fasta_files:
            filename = os.path.basename(filepath)
            print(f"\n  Processing: {filename}")
            try:
                sequences = read_fasta_sequences_only(filepath)
                if not sequences:
                    continue
                text = sequences[0]
                file_res = {'filename': filename, 'sequence_length': len(text), 'patterns': []}
                for pattern in patterns:
                    kmp = KMP(pattern)
                    matches, t = self._measure_time(kmp.search, text)
                    file_res['patterns'].append({'pattern': pattern, 'matches': len(matches), 'time_ms': t*1000, 'first_match_pos': matches[0] if matches else -1})
                    print(f"    Pattern '{pattern}': {len(matches)} matches in {t*1000:.3f} ms")
                results.append(file_res)
            except Exception as e:
                print(f"    Error processing {filename}: {e}")
        return {'benchmark_type': 'real_dataset', 'files_processed': len(results), 'results': results}

    def save_results(self, output_file: str = 'benchmark_results.json') -> None:
        results_with_metadata = {'timestamp': datetime.now().isoformat(), 'benchmarks': self.results}
        with open(output_file, 'w') as f:
            json.dump(results_with_metadata, f, indent=2)
        print(f"\nResults saved to: {output_file}")

    def run_all_benchmarks(self) -> None:
        print("="*70)
        print("KMP ALGORITHM - COMPREHENSIVE BENCHMARK SUITE")
        print("="*70)
        self.results.append(self.test_basic_matching())
        self.results.append(self.test_edge_cases())
        self.results.append(self.benchmark_pattern_length())
        self.results.append(self.benchmark_text_length())
        self.results.append(self.benchmark_multiple_patterns())
        self.results.append(self.benchmark_real_dataset())
        print("\n" + "="*70)
        print("ALL BENCHMARKS COMPLETED")
        print("="*70)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Try to locate dataset at the expected path relative to project root
    dataset_path_candidates = [
        os.path.join(script_dir, '..', 'DnA_dataset', 'ncbi_dataset', 'data'),
        os.path.join(script_dir, '..', '..', 'DnA_dataset', 'ncbi_dataset', 'data')
    ]
    dataset_path = next((p for p in dataset_path_candidates if os.path.exists(p)), None)

    bench = KMPBenchmark(dataset_path=dataset_path)
    bench.run_all_benchmarks()
    output_file = os.path.join(script_dir, 'benchmark_results.json')
    bench.save_results(output_file)

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    total_passed = total_failed = 0
    for result in bench.results:
        if 'passed' in result and 'failed' in result:
            total_passed += result['passed']
            total_failed += result['failed']
    total = total_passed + total_failed
    if total > 0:
        print(f"\nTotal Test Cases: {total}")
        print(f"Passed: {total_passed}")
        print(f"Failed: {total_failed}")
        print(f"Success Rate: {total_passed/total*100:.1f}%")


if __name__ == '__main__':
    main()
