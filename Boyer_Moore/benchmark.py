"""
Benchmarking and Testing Harness for Boyer-Moore Algorithm

This module provides comprehensive testing and benchmarking capabilities
for the Boyer-Moore string matching algorithm on DNA sequences.
"""

import time
import os
import sys
from typing import List, Dict, Tuple
import json
from datetime import datetime

# Add parent directory to path to import our modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from boyer_moore import BoyerMoore, search_multiple_patterns, find_approximate_matches
from utils import (
    read_fasta_file, read_fasta_single_sequence, read_fasta_sequences_only,
    get_all_fasta_files, generate_random_dna, validate_dna_sequence, 
    get_reverse_complement, calculate_gc_content, read_fasta_generator
)


class BoyerMooreBenchmark:
    """
    Comprehensive benchmarking suite for the Boyer-Moore algorithm.
    
    Attributes:
        results (List[Dict]): List of benchmark results
        dataset_path (str): Path to the DNA dataset directory
    """
    
    def __init__(self, dataset_path: str = None):
        """
        Initialize the benchmark suite.
        
        Args:
            dataset_path (str): Path to the DNA dataset directory (optional)
        """
        self.results = []
        self.dataset_path = dataset_path
        
    def _measure_time(self, func, *args, **kwargs) -> Tuple[any, float]:
        """
        Measure the execution time of a function.
        
        Args:
            func: Function to measure
            *args: Positional arguments for the function
            **kwargs: Keyword arguments for the function
        
        Returns:
            Tuple[any, float]: (function result, execution time in seconds)
        """
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        return result, elapsed_time
    
    def test_basic_matching(self) -> Dict:
        """
        Test basic pattern matching functionality.
        
        Returns:
            Dict: Test results with status and details
        """
        print("\n=== Testing Basic Pattern Matching ===")
        test_cases = [
            {
                'name': 'Simple match',
                'text': 'ACGTACGTACGT',
                'pattern': 'ACG',
                'expected': [0, 4, 8]
            },
            {
                'name': 'No match',
                'text': 'ACGTACGTACGT',
                'pattern': 'TTT',
                'expected': []
            },
            {
                'name': 'Single match',
                'text': 'ACGTACGTACGT',
                'pattern': 'GTACG',
                'expected': [2, 6]
            },
            {
                'name': 'Pattern equals text',
                'text': 'ACGT',
                'pattern': 'ACGT',
                'expected': [0]
            },
            {
                'name': 'Overlapping matches',
                'text': 'AAAAAAA',
                'pattern': 'AAA',
                'expected': [0, 1, 2, 3, 4]
            }
        ]
        
        passed = 0
        failed = 0
        details = []
        
        for test in test_cases:
            bm = BoyerMoore(test['pattern'])
            result = bm.search(test['text'])
            
            if result == test['expected']:
                passed += 1
                status = 'PASS'
            else:
                failed += 1
                status = 'FAIL'
            
            details.append({
                'test': test['name'],
                'status': status,
                'expected': test['expected'],
                'got': result
            })
            
            print(f"  {test['name']}: {status}")
        
        print(f"\nBasic Tests: {passed} passed, {failed} failed")
        
        return {
            'test_type': 'basic_matching',
            'passed': passed,
            'failed': failed,
            'details': details
        }
    
    def test_edge_cases(self) -> Dict:
        """
        Test edge cases and boundary conditions.
        
        Returns:
            Dict: Test results with status and details
        """
        print("\n=== Testing Edge Cases ===")
        test_cases = []
        details = []
        passed = 0
        failed = 0
        
        # Test 1: Single character pattern
        try:
            bm = BoyerMoore('A')
            result = bm.search('ACGTACGT')
            expected = [0, 4]
            if result == expected:
                passed += 1
                print("  Single character pattern: PASS")
            else:
                failed += 1
                print("  Single character pattern: FAIL")
        except Exception as e:
            failed += 1
            print(f"  Single character pattern: FAIL ({str(e)})")
        
        # Test 2: Pattern longer than text
        try:
            bm = BoyerMoore('ACGTACGTACGT')
            result = bm.search('ACG')
            expected = []
            if result == expected:
                passed += 1
                print("  Pattern longer than text: PASS")
            else:
                failed += 1
                print("  Pattern longer than text: FAIL")
        except Exception as e:
            failed += 1
            print(f"  Pattern longer than text: FAIL ({str(e)})")
        
        # Test 3: Case insensitivity
        try:
            bm = BoyerMoore('acgt')
            result = bm.search('ACGTACGT')
            expected = [0, 4]
            if result == expected:
                passed += 1
                print("  Case insensitivity: PASS")
            else:
                failed += 1
                print("  Case insensitivity: FAIL")
        except Exception as e:
            failed += 1
            print(f"  Case insensitivity: FAIL ({str(e)})")
        
        # Test 4: Pattern with N (ambiguous base)
        try:
            bm = BoyerMoore('ACNT')
            result = bm.search('ACNTACNT')
            expected = [0, 4]
            if result == expected:
                passed += 1
                print("  Pattern with N: PASS")
            else:
                failed += 1
                print("  Pattern with N: FAIL")
        except Exception as e:
            failed += 1
            print(f"  Pattern with N: FAIL ({str(e)})")
        
        print(f"\nEdge Case Tests: {passed} passed, {failed} failed")
        
        return {
            'test_type': 'edge_cases',
            'passed': passed,
            'failed': failed
        }
    
    def benchmark_pattern_length(self, text_length: int = 100000) -> Dict:
        """
        Benchmark performance with varying pattern lengths.
        
        Args:
            text_length (int): Length of the test text
        
        Returns:
            Dict: Benchmark results
        """
        print("\n=== Benchmarking Pattern Length Impact ===")
        
        # Generate random DNA text
        text = generate_random_dna(text_length, seed=42)
        pattern_lengths = [5, 10, 20, 50, 100, 200]
        
        results = []
        
        for length in pattern_lengths:
            # Create pattern from the text to ensure at least one match
            pattern = text[1000:1000+length]
            
            bm = BoyerMoore(pattern)
            matches, exec_time = self._measure_time(bm.search, text)
            
            results.append({
                'pattern_length': length,
                'text_length': text_length,
                'matches_found': len(matches),
                'time_seconds': exec_time,
                'time_ms': exec_time * 1000
            })
            
            print(f"  Pattern length {length:3d}: {exec_time*1000:8.3f} ms ({len(matches)} matches)")
        
        return {
            'benchmark_type': 'pattern_length',
            'text_length': text_length,
            'results': results
        }
    
    def benchmark_text_length(self, pattern_length: int = 20) -> Dict:
        """
        Benchmark performance with varying text lengths.
        
        Args:
            pattern_length (int): Length of the pattern
        
        Returns:
            Dict: Benchmark results
        """
        print("\n=== Benchmarking Text Length Impact ===")
        
        text_lengths = [10000, 50000, 100000, 500000, 1000000]
        pattern = generate_random_dna(pattern_length, seed=42)
        
        results = []
        
        for length in text_lengths:
            text = generate_random_dna(length, seed=100)
            
            bm = BoyerMoore(pattern)
            matches, exec_time = self._measure_time(bm.search, text)
            
            chars_per_sec = length / exec_time if exec_time > 0 else float('inf')
            
            results.append({
                'text_length': length,
                'pattern_length': pattern_length,
                'matches_found': len(matches),
                'time_seconds': exec_time,
                'time_ms': exec_time * 1000,
                'chars_per_second': chars_per_sec
            })
            
            if exec_time > 0:
                print(f"  Text length {length:7d}: {exec_time*1000:8.3f} ms "
                      f"({length/exec_time/1e6:.2f} M chars/sec)")
            else:
                print(f"  Text length {length:7d}: < 0.001 ms (very fast)")
        
        return {
            'benchmark_type': 'text_length',
            'pattern_length': pattern_length,
            'results': results
        }
    
    def benchmark_real_dataset(self, max_files: int = 5) -> Dict:
        """
        Benchmark on real DNA sequences from the dataset.
        
        Args:
            max_files (int): Maximum number of files to process
        
        Returns:
            Dict: Benchmark results
        """
        print("\n=== Benchmarking on Real DNA Dataset ===")
        
        if not self.dataset_path or not os.path.exists(self.dataset_path):
            print("  Dataset path not found. Skipping real dataset benchmark.")
            return {'benchmark_type': 'real_dataset', 'status': 'skipped'}
        
        # Get FASTA files
        fasta_files = get_all_fasta_files(self.dataset_path, recursive=True)
        
        if not fasta_files:
            print("  No FASTA files found in dataset. Skipping.")
            return {'benchmark_type': 'real_dataset', 'status': 'no_files'}
        
        fasta_files = fasta_files[:max_files]
        print(f"  Processing {len(fasta_files)} files...")
        
        # Common DNA patterns to search for
        patterns = [
            'ATGCATGC',      # 8 bp
            'GCTAGCTA',      # 8 bp
            'TATAAA',        # TATA box
            'CAAT',          # CAAT box
            'GAATTC',        # EcoRI restriction site
            'GGATCC',        # BamHI restriction site
        ]
        
        results = []
        
        for filepath in fasta_files:
            filename = os.path.basename(filepath)
            print(f"\n  Processing: {filename}")
            
            try:
                # Read the sequence
                sequences = read_fasta_sequences_only(filepath)
                if not sequences:
                    continue
                
                text = sequences[0]  # Use first sequence
                text_length = len(text)
                
                file_results = {
                    'filename': filename,
                    'sequence_length': text_length,
                    'patterns': []
                }
                
                for pattern in patterns:
                    bm = BoyerMoore(pattern)
                    matches, exec_time = self._measure_time(bm.search, text)
                    
                    file_results['patterns'].append({
                        'pattern': pattern,
                        'matches': len(matches),
                        'time_ms': exec_time * 1000,
                        'first_match_pos': matches[0] if matches else -1
                    })
                    
                    print(f"    Pattern '{pattern}': {len(matches)} matches in {exec_time*1000:.3f} ms")
                
                results.append(file_results)
                
            except Exception as e:
                print(f"    Error processing {filename}: {str(e)}")
        
        return {
            'benchmark_type': 'real_dataset',
            'files_processed': len(results),
            'results': results
        }
    
    def benchmark_multiple_patterns(self, text_length: int = 100000, num_patterns: int = 10) -> Dict:
        """
        Benchmark searching for multiple patterns in the same text.
        
        Args:
            text_length (int): Length of the test text
            num_patterns (int): Number of patterns to search for
        
        Returns:
            Dict: Benchmark results
        """
        print("\n=== Benchmarking Multiple Pattern Search ===")
        
        text = generate_random_dna(text_length, seed=42)
        
        # Generate patterns of varying lengths
        patterns = []
        for i in range(num_patterns):
            length = 10 + (i * 5)  # Patterns from 10 to 10+5*(num_patterns-1) bp
            pattern = text[i*100:i*100+length]
            patterns.append(pattern)
        
        # Benchmark
        matches_dict, exec_time = self._measure_time(search_multiple_patterns, text, patterns)
        
        total_matches = sum(len(matches) for matches in matches_dict.values())
        
        print(f"  Searched for {num_patterns} patterns in text of length {text_length}")
        print(f"  Total time: {exec_time*1000:.3f} ms")
        print(f"  Total matches found: {total_matches}")
        print(f"  Average time per pattern: {exec_time*1000/num_patterns:.3f} ms")
        
        pattern_results = []
        for pattern, matches in matches_dict.items():
            pattern_results.append({
                'pattern': pattern[:20] + '...' if len(pattern) > 20 else pattern,
                'pattern_length': len(pattern),
                'matches': len(matches)
            })
        
        return {
            'benchmark_type': 'multiple_patterns',
            'text_length': text_length,
            'num_patterns': num_patterns,
            'total_matches': total_matches,
            'total_time_ms': exec_time * 1000,
            'avg_time_per_pattern_ms': exec_time * 1000 / num_patterns,
            'pattern_results': pattern_results
        }
    
    def save_results(self, output_file: str = 'benchmark_results.json') -> None:
        """
        Save benchmark results to a JSON file.
        
        Args:
            output_file (str): Output file path
        """
        results_with_metadata = {
            'timestamp': datetime.now().isoformat(),
            'benchmarks': self.results
        }
        
        with open(output_file, 'w') as f:
            json.dump(results_with_metadata, f, indent=2)
        
        print(f"\nResults saved to: {output_file}")
    
    def run_all_benchmarks(self) -> None:
        """
        Run all benchmarks and tests.
        """
        print("="*70)
        print("BOYER-MOORE ALGORITHM - COMPREHENSIVE BENCHMARK SUITE")
        print("="*70)
        
        # Tests
        self.results.append(self.test_basic_matching())
        self.results.append(self.test_edge_cases())
        
        # Benchmarks
        self.results.append(self.benchmark_pattern_length())
        self.results.append(self.benchmark_text_length())
        self.results.append(self.benchmark_multiple_patterns())
        
        # Real dataset (if available)
        self.results.append(self.benchmark_real_dataset())
        
        print("\n" + "="*70)
        print("ALL BENCHMARKS COMPLETED")
        print("="*70)


def main():
    """
    Main function to run benchmarks.
    """
    # Get the dataset path (relative to this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_path = os.path.join(script_dir, '..', 'STARK', 'Suffix Arrays-Trees', 
                                 'DnA_dataset', 'ncbi_dataset', 'data')
    
    # Create benchmark suite
    benchmark = BoyerMooreBenchmark(dataset_path=dataset_path)
    
    # Run all benchmarks
    benchmark.run_all_benchmarks()
    
    # Save results
    output_file = os.path.join(script_dir, 'benchmark_results.json')
    benchmark.save_results(output_file)
    
    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    total_tests = 0
    total_passed = 0
    total_failed = 0
    
    for result in benchmark.results:
        if 'passed' in result and 'failed' in result:
            total_tests += 1
            total_passed += result['passed']
            total_failed += result['failed']
    
    print(f"\nTotal Test Cases: {total_passed + total_failed}")
    print(f"Passed: {total_passed}")
    print(f"Failed: {total_failed}")
    print(f"Success Rate: {total_passed/(total_passed+total_failed)*100:.1f}%" 
          if (total_passed + total_failed) > 0 else "N/A")


if __name__ == '__main__':
    main()
