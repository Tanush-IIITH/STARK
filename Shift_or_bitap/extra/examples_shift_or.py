# examples_shift_or.py
from shift_or_bitap import ShiftOr

def smoke_tests():
    text = "ACGTACGTATGCATGCACGT"
    pattern = "ATGC"
    so = ShiftOr(pattern)
    print("exact matches:", so.search(text))
    so_approx = ShiftOr("ATGG", max_errors=1)
    print("approx matches (k=1):", so_approx.search_approx(text, 1))

if __name__ == "__main__":
    smoke_tests()
