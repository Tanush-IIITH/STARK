"""
Utility Functions for DNA Sequence Processing (Shift-Or/Bitap)

This module duplicates the utilities used by the Boyer-Moore module for
reading and processing DNA sequences, to keep the KMP module self-contained.

Author: DNA Pattern Matching Project
Date: November 2025
"""

import os
from typing import List, Tuple, Dict, Generator


def read_fasta_file(filepath: str) -> Dict[str, str]:
    """
    Read a FASTA file and return sequences with their headers.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    sequences: Dict[str, str] = {}
    current_header = None
    current_sequence: List[str] = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_sequence)
                current_header = line[1:]
                current_sequence = []
            else:
                if current_header is None:
                    raise ValueError("Sequence data found before header in FASTA file")
                current_sequence.append(line.upper())
        if current_header is not None:
            sequences[current_header] = ''.join(current_sequence)
    return sequences


def read_fasta_sequences_only(filepath: str) -> List[str]:
    sequences_dict = read_fasta_file(filepath)
    return list(sequences_dict.values())


def read_fasta_single_sequence(filepath: str) -> str:
    sequences = read_fasta_sequences_only(filepath)
    return ''.join(sequences)


def read_fasta_generator(filepath: str) -> Generator[Tuple[str, str], None, None]:
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    current_header = None
    current_sequence: List[str] = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    yield (current_header, ''.join(current_sequence))
                current_header = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line.upper())
        if current_header is not None:
            yield (current_header, ''.join(current_sequence))


def get_all_fasta_files(directory: str, recursive: bool = True) -> List[str]:
    fasta_extensions = ['.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn']
    fasta_files: List[str] = []

    if recursive:
        for root, _dirs, files in os.walk(directory):
            for file in files:
                if any(file.lower().endswith(ext) for ext in fasta_extensions):
                    fasta_files.append(os.path.join(root, file))
    else:
        for file in os.listdir(directory):
            filepath = os.path.join(directory, file)
            if os.path.isfile(filepath) and any(file.lower().endswith(ext) for ext in fasta_extensions):
                fasta_files.append(filepath)

    return sorted(fasta_files)


def validate_dna_sequence(sequence: str) -> bool:
    valid_chars = set('ACGTN')
    return all(c in valid_chars for c in sequence.upper())


def calculate_gc_content(sequence: str) -> float:
    sequence = sequence.upper()
    if not sequence:
        return 0.0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100


def get_reverse_complement(sequence: str) -> str:
    complement_map = {
        'A': 'T', 'T': 'A',
        'G': 'C', 'C': 'G',
        'N': 'N'
    }
    sequence = sequence.upper()
    complement = ''.join(complement_map.get(base, 'N') for base in sequence)
    return complement[::-1]


def generate_random_dna(length: int, seed: int = None) -> str:
    import random
    if seed is not None:
        random.seed(seed)
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))


def extract_subsequence(sequence: str, start: int, end: int) -> str:
    return sequence[start:end]


def count_nucleotides(sequence: str) -> Dict[str, int]:
    sequence = sequence.upper()
    return {
        'A': sequence.count('A'),
        'C': sequence.count('C'),
        'G': sequence.count('G'),
        'T': sequence.count('T'),
        'N': sequence.count('N')
    }


def split_sequence_into_chunks(sequence: str, chunk_size: int, overlap: int = 0) -> List[str]:
    chunks: List[str] = []
    step = max(1, chunk_size - overlap)
    for i in range(0, len(sequence), step):
        chunk = sequence[i:i + chunk_size]
        if len(chunk) >= chunk_size:
            chunks.append(chunk)
        elif i == 0 or len(chunk) > overlap:
            chunks.append(chunk)
    return chunks


def format_sequence_pretty(sequence: str, line_length: int = 80) -> str:
    lines: List[str] = []
    for i in range(0, len(sequence), line_length):
        lines.append(sequence[i:i + line_length])
    return '\n'.join(lines)
