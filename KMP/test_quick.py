"""
Quick unit-like checks for KMP implementation.

Run: python test_quick.py
"""

from kmp import KMP


def run_tests():
    tests = [
        ("ACGTACGTACGT", "ACG", [0, 4, 8]),
        ("ACGTACGTACGT", "TTT", []),
        ("AAAAAAA", "AAA", [0, 1, 2, 3, 4]),
        ("ACGT", "ACGT", [0]),
        ("", "A", []),
    ]

    passed = 0
    failed = 0

    for text, pattern, expected in tests:
        try:
            if not pattern:
                print(f"Skipping empty pattern test: text={text!r}")
                continue
            k = KMP(pattern)
            res = k.search(text)
            if res == expected:
                print(f"PASS: pattern={pattern!r} text={text!r} -> {res}")
                passed += 1
            else:
                print(f"FAIL: pattern={pattern!r} text={text!r} -> got {res}, expected {expected}")
                failed += 1
        except Exception as e:
            print(f"ERROR: pattern={pattern!r} text={text!r} -> {e}")
            failed += 1

    print(f"\nQuick tests: {passed} passed, {failed} failed")


if __name__ == '__main__':
    run_tests()
