#!/usr/bin/env python3
"""
Generate a synthetic FASTA dataset and inject known pattern occurrences for
testing string-search algorithms (KMP). Produces:

 - dataset/<PREFIX>_genomic.fna
 - dataset/<PREFIX>_patterns.tsv   (tab-separated: pattern<TAB>pos1,pos2,...)

Usage examples:
  python synthetic_kmp.py --length 100000 --pattern ATGCGT --pattern-count 50 --seed 42
  python synthetic_kmp.py --length 50000 --patterns-file my_patterns.txt

Patterns may be provided via --pattern option (repeatable) or --patterns-file
(one pattern per line).
"""
from __future__ import annotations

import argparse
import os
import random
import string
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Tuple, Dict

NUCLEOTIDES = "ACGT"
DEFAULT_WRAP = 70
DEFAULT_LENGTH = 200_000


def wrap_lines(sequence: str, width: int) -> Iterable[str]:
    for idx in range(0, len(sequence), width):
        yield sequence[idx: idx + width]


def generate_sequence(length: int, seed: int | None = None) -> List[str]:
    rng = random.Random(seed)
    return [rng.choice(NUCLEOTIDES) for _ in range(length)]


def inject_patterns(
    seq: List[str],
    patterns: List[str],
    counts: Dict[str, int],
    rng: random.Random,
    allow_overlap: bool,
    min_distance: int,
) -> Dict[str, List[int]]:
    """
    Inject patterns into seq (list of chars). Return mapping pattern -> list of 1-based positions.
    If allow_overlap is False, respect min_distance between insertions.
    """
    n = len(seq)
    occupied = [False] * n  # mark bases that are already used (if no-overlap)
    positions: Dict[str, List[int]] = {p: [] for p in patterns}

    for p in patterns:
        p_len = len(p)
        attempts = 0
        target = counts.get(p, 0)
        while len(positions[p]) < target:
            attempts += 1
            if attempts > target * 1000:
                # too many failures: relax constraints (allow overlap) rather than loop forever
                if not allow_overlap:
                    # allow overlap as fallback
                    allow_overlap = True
                else:
                    raise RuntimeError(f"Could not place pattern {p} {target} times (tried {attempts})")
            pos = rng.randint(0, n - p_len)
            # check occupancy constraints
            if not allow_overlap:
                ok = True
                # ensure min_distance around position
                start = max(0, pos - min_distance)
                end = min(n, pos + p_len + min_distance)
                if any(occupied[i] for i in range(start, end)):
                    ok = False
                if not ok:
                    continue
                # mark occupied
                for i in range(pos, pos + p_len):
                    occupied[i] = True

            # inject the pattern
            seq[pos: pos + p_len] = list(p)
            positions[p].append(pos + 1)  # 1-based positions for humans
    # keep positions sorted
    for p in positions:
        positions[p].sort()
    return positions


def parse_patterns_from_file(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8") as fh:
        pats = [line.strip().upper() for line in fh if line.strip()]
    return pats


def main() -> None:
    parser = argparse.ArgumentParser(description="Create synthetic FASTA + pattern ground truth for KMP tests")
    parser.add_argument("--length", type=int, default=DEFAULT_LENGTH, help="total bases to generate")
    parser.add_argument("--wrap", type=int, default=DEFAULT_WRAP, help="FASTA wrap width")
    parser.add_argument("--seed", type=int, default=None, help="RNG seed")
    parser.add_argument("--prefix", default=None, help="dataset filename prefix")
    parser.add_argument("--pattern", action="append", default=[], help="Pattern to inject (can repeat)")
    parser.add_argument("--pattern-count", type=int, default=10, help="Default count for each supplied pattern")
    parser.add_argument("--patterns-file", type=str, help="File with one pattern per line")
    parser.add_argument("--allow-overlap", action="store_true", help="Allow overlapping pattern placements")
    parser.add_argument("--min-distance", type=int, default=5, help="Min distance between injected pattern occurrences (when overlap not allowed)")
    args = parser.parse_args()

    if args.length <= 0:
        raise ValueError("length must be positive")
    if args.wrap <= 0:
        raise ValueError("wrap must be positive")

    rng = random.Random(args.seed)

    # gather patterns
    patterns = list(args.pattern)
    if args.patterns_file:
        pfromfile = parse_patterns_from_file(Path(args.patterns_file))
        patterns.extend(pfromfile)
    if not patterns:
        # default pattern set (short motifs)
        patterns = ["ATG", "TATA", "GCGC", "AATT", "CCGGA"]

    # build counts dict (same count for every pattern unless user supplied duplicates)
    counts = {}
    for p in patterns:
        counts[p] = counts.get(p, 0) + args.pattern_count

    # prepare sequence
    seq_list = generate_sequence(args.length, seed=args.seed)

    # inject patterns (this mutates seq_list)
    positions = inject_patterns(
        seq_list,
        patterns=list(counts.keys()),
        counts=counts,
        rng=rng,
        allow_overlap=args.allow_overlap,
        min_distance=args.min_distance,
    )

    # prefix/safe
    if args.prefix:
        safe_prefix = "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in args.prefix)
    else:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        safe_prefix = f"SYNTHKMP_{stamp}_{''.join(rng.choices(string.ascii_uppercase, k=4))}"

    script_dir = Path(__file__).resolve().parent
    dataset_dir = script_dir / "dataset"
    dataset_dir.mkdir(parents=True, exist_ok=True)

    fasta_name = dataset_dir / f"{safe_prefix}_genomic.fna"
    patterns_name = dataset_dir / f"{safe_prefix}_patterns.tsv"

    if fasta_name.exists() or patterns_name.exists():
        raise FileExistsError(f"Refusing to overwrite existing dataset: {fasta_name} or {patterns_name}")

    # write FASTA
    header = f">{safe_prefix} synthetic genome (length={len(seq_list)})"
    with fasta_name.open("w", encoding="utf-8") as fh:
        fh.write(header + "\n")
        for line in wrap_lines("".join(seq_list), args.wrap):
            fh.write(line + "\n")

    # write patterns TSV: pattern<TAB>comma-separated positions
    with patterns_name.open("w", encoding="utf-8") as fh:
        fh.write("pattern\tpositions\n")
        for p, poslist in positions.items():
            fh.write(f"{p}\t{','.join(str(x) for x in poslist)}\n")

    print(f"Wrote FASTA: {fasta_name}")
    print(f"Wrote pattern truth: {patterns_name}")
    print("Patterns counts:")
    for p, poslist in positions.items():
        print(f"  {p}: {len(poslist)} occurrences")

if __name__ == "__main__":
    main()
