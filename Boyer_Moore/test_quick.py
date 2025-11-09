"""
Quick Test Script for Boyer-Moore Algorithm

A simple script to quickly verify the implementation is working correctly.

Author: DNA Pattern Matching Project
Date: November 2025
"""

import sys
import os

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from boyer_moore import BoyerMoore


def run_quick_tests():
    """Run quick validation tests."""
    print("\n" + "="*60)
    print("BOYER-MOORE ALGORITHM - QUICK VALIDATION TEST")
    print("="*60 + "\n")
    
    tests_passed = 0
    tests_failed = 0
    
    # Test 1: Basic match
    print("Test 1: Basic pattern matching...")
    bm = BoyerMoore("ATGC")
    result = bm.search("ACGTACGTATGCATGCACGT")
    expected = [8, 12]
    if result == expected:
        print("  ✓ PASS")
        tests_passed += 1
    else:
        print(f"  ✗ FAIL - Expected {expected}, got {result}")
        tests_failed += 1
    
    # Test 2: No match
    print("Test 2: No match scenario...")
    bm = BoyerMoore("ZZZZZ")
    result = bm.search("ACGTACGTATGCATGCACGT")
    expected = []
    if result == expected:
        print("  ✓ PASS")
        tests_passed += 1
    else:
        print(f"  ✗ FAIL - Expected {expected}, got {result}")
        tests_failed += 1
    
    # Test 3: Overlapping matches
    print("Test 3: Overlapping matches...")
    bm = BoyerMoore("AAA")
    result = bm.search("AAAAAAA")
    expected = [0, 1, 2, 3, 4]
    if result == expected:
        print("  ✓ PASS")
        tests_passed += 1
    else:
        print(f"  ✗ FAIL - Expected {expected}, got {result}")
        tests_failed += 1
    
    # Test 4: Case insensitivity
    print("Test 4: Case insensitivity...")
    bm = BoyerMoore("atgc")
    result = bm.search("ATGCATGC")
    expected = [0, 4]
    if result == expected:
        print("  ✓ PASS")
        tests_passed += 1
    else:
        print(f"  ✗ FAIL - Expected {expected}, got {result}")
        tests_failed += 1
    
    # Test 5: First occurrence
    print("Test 5: Finding first occurrence...")
    bm = BoyerMoore("ATGC")
    result = bm.search_first("ACGTATGCATGCACGT")
    expected = 4
    if result == expected:
        print("  ✓ PASS")
        tests_passed += 1
    else:
        print(f"  ✗ FAIL - Expected {expected}, got {result}")
        tests_failed += 1
    
    # Test 6: Count matches
    print("Test 6: Counting matches...")
    bm = BoyerMoore("AT")
    result = bm.count_matches("ATATATATAT")
    expected = 5
    if result == expected:
        print("  ✓ PASS")
        tests_passed += 1
    else:
        print(f"  ✗ FAIL - Expected {expected}, got {result}")
        tests_failed += 1
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    print(f"Tests Passed: {tests_passed}")
    print(f"Tests Failed: {tests_failed}")
    print(f"Success Rate: {tests_passed/(tests_passed+tests_failed)*100:.1f}%")
    
    if tests_failed == 0:
        print("\n✓ All tests passed! Implementation is working correctly.")
        return 0
    else:
        print(f"\n✗ {tests_failed} test(s) failed. Please check the implementation.")
        return 1


if __name__ == '__main__':
    exit_code = run_quick_tests()
    sys.exit(exit_code)
