"""Generate a synthetic FASTA dataset that mirrors the structure of the
existing genomes under ``dataset/``.

Each run emits a new ``*.fna`` file whose header and line wrapping follow the
same conventions as NCBI genomes bundled with the project.
"""
from __future__ import annotations

import argparse
import os
import random
import string
from datetime import datetime, UTC
from pathlib import Path
from typing import Iterable

NUCLEOTIDES = "ACGT"
DEFAULT_WRAP = 70
DEFAULT_LENGTH = 200_000


def wrap_lines(sequence: str, width: int) -> Iterable[str]:
    """Yield ``sequence`` broken into ``width``-sized chunks."""
    for idx in range(0, len(sequence), width):
        yield sequence[idx : idx + width]


def generate_sequence(length: int, seed: int | None = None) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(NUCLEOTIDES) for _ in range(length))


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a synthetic FASTA dataset")
    parser.add_argument(
        "--length",
        type=int,
        default=DEFAULT_LENGTH,
        help="Number of bases to generate (default: %(default)s)",
    )
    parser.add_argument(
        "--wrap",
        type=int,
        default=DEFAULT_WRAP,
        help="Line width for FASTA output (default: %(default)s)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional RNG seed for reproducible output",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Custom dataset prefix (omit suffix). Defaults to timestamped SYNTH name.",
    )
    args = parser.parse_args()

    if args.length <= 0:
        raise ValueError("length must be positive")
    if args.wrap <= 0:
        raise ValueError("wrap must be positive")

    script_dir = Path(__file__).resolve().parent
    dataset_dir = script_dir / "dataset"
    dataset_dir.mkdir(parents=True, exist_ok=True)

    if args.prefix:
        safe_prefix = "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in args.prefix)
    else:
        stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S")
        safe_prefix = f"SYNTH_{stamp}_{''.join(random.choices(string.ascii_uppercase, k=4))}"

    filename = f"{safe_prefix}_genomic.fna"
    output_path = dataset_dir / filename

    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing dataset: {output_path}")

    sequence = generate_sequence(args.length, seed=args.seed)
    header = f">{safe_prefix} synthetic genome (length={len(sequence)})"

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(header + "\n")
        for line in wrap_lines(sequence, args.wrap):
            handle.write(line + "\n")

    print(f"Synthetic dataset written to {output_path}")


if __name__ == "__main__":
    main()
