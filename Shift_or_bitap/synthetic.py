"""
Generate synthetic FASTA datasets for DNA pattern matching benchmarks.

This script generates random DNA sequences in FASTA format with customizable
parameters for benchmarking pattern matching algorithms.

Compatible with Python 3.10+
"""

import argparse
import os
import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

NUCLEOTIDES = "ACGT"
DEFAULT_WRAP = 70
DEFAULT_LENGTH = 200_000


def wrap_lines(sequence: str, width: int) -> Iterable[str]:
    """Yield sequence broken into width-sized chunks."""
    for idx in range(0, len(sequence), width):
        yield sequence[idx : idx + width]


def generate_sequence(length: int, seed: int = None) -> str:
    """
    Generate a random DNA sequence of given length.

    Args:
        length: Number of base pairs to generate
        seed: Random seed for reproducibility (optional)

    Returns:
        String of random DNA nucleotides (A, C, G, T)
    """
    rng = random.Random(seed)
    return "".join(rng.choice(NUCLEOTIDES) for _ in range(length))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate synthetic FASTA dataset for DNA pattern matching",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate 100K bp dataset
  python synthetic.py --length 100000 --seed 42 --name SYNTH_100K \
      --output synthetic_datasets/SYNTH_100K/SYNTH_100K.fna

  # Generate 1M bp dataset with custom prefix
  python synthetic.py --length 1000000 --seed 46 --prefix SYNTH_1M

  # Use default settings
  python synthetic.py --prefix TEST
        """
    )

    parser.add_argument(
        "--length",
        type=int,
        default=DEFAULT_LENGTH,
        help=f"Length of the synthetic sequence in base pairs (default: {DEFAULT_LENGTH:,})",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility (default: uses system random)",
    )

    parser.add_argument(
        "--name",
        type=str,
        default=None,
        help="Name for the dataset (used in FASTA header and default output path)",
    )

    parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="Alternative to --name (for backward compatibility)",
    )

    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output FASTA file path (default: synthetic_datasets/<name>/<name>.fna)",
    )

    parser.add_argument(
        "--wrap",
        type=int,
        default=DEFAULT_WRAP,
        help=f"Number of characters per line in FASTA (default: {DEFAULT_WRAP})",
    )

    args = parser.parse_args()

    # Determine dataset name (prefer --name, fall back to --prefix)
    dataset_name = args.name or args.prefix
    if not dataset_name:
        dataset_name = "SYNTHETIC"

    # Determine output path
    if args.output:
        output_path = Path(args.output)
    else:
        # Default: synthetic_datasets/<name>/<name>.fna
        output_path = Path("synthetic_datasets") / dataset_name / f"{dataset_name}.fna"

    # Create output directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Generate sequence
    print(f"Generating {args.length:,} bp synthetic DNA sequence...")
    print(f"  Dataset name: {dataset_name}")
    print(f"  Random seed: {args.seed if args.seed is not None else 'random'}")

    sequence = generate_sequence(args.length, args.seed)

    # Write FASTA file
    print(f"Writing to {output_path}...")

    # Use timezone.utc instead of UTC for Python 3.10 compatibility
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    with open(output_path, "w") as f:
        # FASTA header
        seed_str = f"seed={args.seed}" if args.seed is not None else "seed=random"
        header = f">{dataset_name} length={args.length} {seed_str} generated={timestamp}"
        f.write(header + "\n")

        # Write sequence with line wrapping
        for line in wrap_lines(sequence, args.wrap):
            f.write(line + "\n")

    # Print summary
    file_size_kb = output_path.stat().st_size / 1024
    print(f"\n✓ Successfully generated synthetic dataset!")
    print(f"  Output file: {output_path}")
    print(f"  Sequence length: {args.length:,} bp")
    print(f"  File size: {file_size_kb:.2f} KB")

    # Calculate composition
    a_count = sequence.count('A')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    t_count = sequence.count('T')
    gc_content = (g_count + c_count) / args.length * 100

    print(f"\n  Nucleotide composition:")
    print(f"    A: {a_count:,} ({a_count/args.length*100:.1f}%)")
    print(f"    C: {c_count:,} ({c_count/args.length*100:.1f}%)")
    print(f"    G: {g_count:,} ({g_count/args.length*100:.1f}%)")
    print(f"    T: {t_count:,} ({t_count/args.length*100:.1f}%)")
    print(f"    GC content: {gc_content:.1f}%")


if __name__ == "__main__":
    main()
