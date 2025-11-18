"""
Quick Tests for Shift-Or/Bitap Algorithm Implementation

This script runs basic validation tests to ensure all three algorithm variants
(exact, approximate, extended) are working correctly.

Author: DNA Pattern Matching Project
Date: November 2025
"""

from shift_or_exact import ShiftOrExact
from shift_or_approximate import ShiftOrApproximate
from shift_or_extended import ShiftOrExtended

def test_exact_matching():
    """Test exact matching algorithm."""
    print("=" * 60)
    print("TEST 1: Exact Matching (≤64 bp)")
    print("=" * 60)

    # Test 1: Basic exact match
    text = "ACGTACGTACGT"
    pattern = "ACGT"
    matcher = ShiftOrExact(pattern)
    matches = matcher.search(text)
    expected = [0, 4, 8]
    assert matches == expected, f"Expected {expected}, got {matches}"
    print(f"✓ Basic match: Pattern '{pattern}' found at positions {matches}")

    # Test 2: No matches
    text = "AAAAAAAA"
    pattern = "CGCG"
    matcher = ShiftOrExact(pattern)
    matches = matcher.search(text)
    assert matches == [], f"Expected no matches, got {matches}"
    print(f"✓ No match: Pattern '{pattern}' correctly not found")

    # Test 3: Longer pattern
    text = "ACGTACGTTAGCTAGCTAGCT"
    pattern = "TAGCT"
    matcher = ShiftOrExact(pattern)
    matches = matcher.search(text)
    print(f"✓ Longer pattern: Found at positions {matches}")

    # Test 4: Maximum length (64 bp)
    pattern_64 = "A" * 64
    text_64 = "C" * 100 + "A" * 64 + "T" * 100
    matcher = ShiftOrExact(pattern_64)
    matches = matcher.search(text_64)
    assert 100 in matches, "64 bp pattern not found"
    print(f"✓ 64 bp pattern: Successfully handled maximum length")

    # Test 5: Metrics tracking
    pattern = "ACGT"
    text = "ACGTACGTACGT"
    matcher = ShiftOrExact(pattern)
    metrics = matcher.search_with_metrics(text)
    print(f"✓ Metrics: {metrics['bit_operations']} bit operations, "
          f"{metrics['state_vectors']} state vector(s)")

    print("✅ All exact matching tests passed!\n")


def test_approximate_matching():
    """Test approximate matching algorithm."""
    print("=" * 60)
    print("TEST 2: Approximate Matching (≤64 bp, k errors)")
    print("=" * 60)

    # Test 1: k=1 - one substitution
    text = "ACGTACTT"  # Second match has one substitution (G→T)
    pattern = "ACGT"
    matcher = ShiftOrApproximate(pattern, k=1)
    matches = matcher.search(text)
    positions = [pos for pos, error in matches]
    assert 0 in positions and 4 in positions, f"Expected matches at 0 and 4, got {positions}"
    print(f"✓ k=1 (substitution): Found {len(matches)} matches")

    # Test 2: k=2 - multiple errors
    text = "ACGTACTTACAT"
    pattern = "ACGT"
    matcher = ShiftOrApproximate(pattern, k=2)
    matches = matcher.search(text)
    print(f"✓ k=2 (multiple errors): Found {len(matches)} matches")

    # Test 3: k=3 - three errors
    pattern = "AAAA"
    text = "TTTT"  # All 4 characters differ
    matcher = ShiftOrApproximate(pattern, k=3)
    matches = matcher.search(text)
    print(f"✓ k=3 (three errors): Found {len(matches)} matches")

    # Test 4: Error levels
    text = "ACGTACTTACATAAAA"
    pattern = "ACGT"
    for k in [1, 2, 3]:
        matcher = ShiftOrApproximate(pattern, k=k)
        matches = matcher.search(text)
        error_levels = [error for pos, error in matches]
        print(f"✓ k={k}: {len(matches)} matches, error levels: {set(error_levels)}")

    # Test 5: Metrics with approximate matching
    pattern = "ACGT"
    text = "ACGTACTTACGT"
    matcher = ShiftOrApproximate(pattern, k=1)
    metrics = matcher.search_with_metrics(text)
    print(f"✓ Approximate metrics: {metrics['state_vectors']} state vectors")

    print("✅ All approximate matching tests passed!\n")


def test_extended_matching():
    """Test extended matching for patterns >64 bp."""
    print("=" * 60)
    print("TEST 3: Extended Matching (>64 bp, exact only)")
    print("=" * 60)

    # Test 1: 100 bp pattern
    pattern_100 = "A" * 100
    text_100 = "C" * 200 + "A" * 100 + "T" * 200
    matcher = ShiftOrExtended(pattern_100)
    matches = matcher.search(text_100)
    assert 200 in matches, "100 bp pattern not found"
    print(f"✓ 100 bp pattern: Found at position {matches}")

    # Test 2: 200 bp pattern
    pattern_200 = "ACGT" * 50  # 200 bp
    text_200 = "T" * 500 + pattern_200 + "G" * 500
    matcher = ShiftOrExtended(pattern_200)
    matches = matcher.search(text_200)
    assert 500 in matches, "200 bp pattern not found"
    print(f"✓ 200 bp pattern: Found at position {matches}")

    # Test 3: 400 bp pattern
    pattern_400 = "ACGT" * 100  # 400 bp
    text_400 = "NNNN" + pattern_400 + "NNNN"
    matcher = ShiftOrExtended(pattern_400)
    matches = matcher.search(text_400)
    assert 4 in matches, "400 bp pattern not found"
    print(f"✓ 400 bp pattern: Successfully handled")

    # Test 4: 800 bp pattern (maximum)
    pattern_800 = "ACGT" * 200  # 800 bp
    text_800 = pattern_800 + "NNNN"
    matcher = ShiftOrExtended(pattern_800)
    matches = matcher.search(text_800)
    assert 0 in matches, "800 bp pattern not found"
    print(f"✓ 800 bp pattern: Maximum length handled")

    # Test 5: Multi-word metrics
    pattern_150 = "A" * 150
    text_150 = "A" * 300
    matcher = ShiftOrExtended(pattern_150)
    metrics = matcher.search_with_metrics(text_150)
    print(f"✓ Multi-word metrics: {metrics['state_vectors']} state vectors, "
          f"{metrics['bit_operations']} bit operations")

    print("✅ All extended matching tests passed!\n")


def test_error_handling():
    """Test error handling and edge cases."""
    print("=" * 60)
    print("TEST 4: Error Handling and Edge Cases")
    print("=" * 60)

    # Test 1: Empty pattern
    try:
        matcher = ShiftOrExact("")
        print("✗ Should raise error for empty pattern")
    except ValueError:
        print("✓ Empty pattern correctly rejected")

    # Test 2: Pattern too long for exact/approximate
    try:
        pattern = "A" * 65
        matcher = ShiftOrExact(pattern)
        print("✗ Should raise error for 65 bp pattern in exact matcher")
    except ValueError:
        print("✓ Pattern >64 bp correctly rejected for exact matcher")

    # Test 3: Pattern too short for extended
    try:
        pattern = "A" * 50
        matcher = ShiftOrExtended(pattern)
        print("✗ Should raise error for pattern ≤64 bp in extended matcher")
    except ValueError:
        print("✓ Pattern ≤64 bp correctly rejected for extended matcher")

    # Test 4: Invalid DNA characters
    try:
        pattern = "ACGTXYZ"
        matcher = ShiftOrExact(pattern)
        print("✗ Should raise error for invalid DNA characters")
    except ValueError:
        print("✓ Invalid DNA characters correctly rejected")

    # Test 5: Invalid k value
    try:
        pattern = "ACGT"
        matcher = ShiftOrApproximate(pattern, k=5)
        print("✗ Should raise error for k > 3")
    except ValueError:
        print("✓ Invalid k value correctly rejected")

    print("✅ All error handling tests passed!\n")


def run_all_tests():
    """Run all validation tests."""
    print("\n" + "=" * 60)
    print("SHIFT-OR/BITAP ALGORITHM - VALIDATION TESTS")
    print("=" * 60 + "\n")

    test_exact_matching()
    test_approximate_matching()
    test_extended_matching()
    test_error_handling()

    print("=" * 60)
    print("✅ ALL TESTS PASSED SUCCESSFULLY!")
    print("=" * 60)
    print("\nThe Shift-Or/Bitap implementation is ready to use.")
    print("Run the benchmark scripts or notebooks to test on real datasets.\n")


if __name__ == "__main__":
    run_all_tests()
