"""
Utility Functions for DNA Sequence Processing

This module provides helper functions for reading and processing DNA sequences
from FASTA files and other formats commonly used in bioinformatics.
"""

import os
from typing import List, Tuple, Dict, Generator
import re


def read_fasta_file(filepath: str) -> Dict[str, str]:
    """
    Read a FASTA file and return sequences with their headers.
    
    FASTA format consists of:
    - Header lines starting with '>'
    - Sequence lines following each header
    
    Args:
        filepath (str): Path to the FASTA file
    
    Returns:
        Dict[str, str]: Dictionary mapping sequence headers to their sequences
    
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    sequences = {}
    current_header = None
    current_sequence = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            if not line:  # Skip empty lines
                continue
            
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None:
                    sequences[current_header] = ''.join(current_sequence)
                
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_sequence = []
            else:
                # Add to current sequence
                if current_header is None:
                    raise ValueError("Sequence data found before header in FASTA file")
                current_sequence.append(line.upper())
        
        # Save the last sequence
        if current_header is not None:
            sequences[current_header] = ''.join(current_sequence)
    
    return sequences


def read_fasta_sequences_only(filepath: str) -> List[str]:
    """
    Read a FASTA file and return only the sequences (without headers).
    
    Args:
        filepath (str): Path to the FASTA file
    
    Returns:
        List[str]: List of DNA sequences
    """
    sequences_dict = read_fasta_file(filepath)
    return list(sequences_dict.values())


def read_fasta_single_sequence(filepath: str) -> str:
    """
    Read a FASTA file and concatenate all sequences into one string.
    
    Useful for files with a single sequence or when treating multiple
    sequences as one continuous sequence.
    
    Args:
        filepath (str): Path to the FASTA file
    
    Returns:
        str: Concatenated DNA sequence
    """
    sequences = read_fasta_sequences_only(filepath)
    return ''.join(sequences)


def read_fasta_generator(filepath: str) -> Generator[Tuple[str, str], None, None]:
    """
    Read a FASTA file lazily using a generator for memory efficiency.
    
    Useful for processing very large FASTA files without loading
    the entire file into memory.
    
    Args:
        filepath (str): Path to the FASTA file
    
    Yields:
        Tuple[str, str]: Tuple of (header, sequence)
    
    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    current_header = None
    current_sequence = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            if not line:
                continue
            
            if line.startswith('>'):
                # Yield previous sequence if exists
                if current_header is not None:
                    yield (current_header, ''.join(current_sequence))
                
                # Start new sequence
                current_header = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line.upper())
        
        # Yield the last sequence
        if current_header is not None:
            yield (current_header, ''.join(current_sequence))


def get_all_fasta_files(directory: str, recursive: bool = True) -> List[str]:
    """
    Get all FASTA files in a directory.
    
    Args:
        directory (str): Path to the directory
        recursive (bool): Whether to search recursively in subdirectories
    
    Returns:
        List[str]: List of paths to FASTA files
    """
    fasta_extensions = ['.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn']
    fasta_files = []
    
    if recursive:
        for root, dirs, files in os.walk(directory):
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
    """
    Validate if a string is a valid DNA sequence.
    
    Valid DNA sequences contain only A, C, G, T, and N characters.
    N represents any nucleotide (used for ambiguous bases).
    
    Args:
        sequence (str): The sequence to validate
    
    Returns:
        bool: True if valid DNA sequence, False otherwise
    """
    valid_chars = set('ACGTN')
    return all(c in valid_chars for c in sequence.upper())


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate the GC content of a DNA sequence.
    
    GC content is the percentage of bases that are G or C.
    
    Args:
        sequence (str): DNA sequence
    
    Returns:
        float: GC content as a percentage (0-100)
    """
    sequence = sequence.upper()
    if not sequence:
        return 0.0
    
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100


def get_reverse_complement(sequence: str) -> str:
    """
    Get the reverse complement of a DNA sequence.
    
    The complement pairs are: A-T and G-C
    
    Args:
        sequence (str): DNA sequence
    
    Returns:
        str: Reverse complement of the sequence
    """
    complement_map = {
        'A': 'T', 'T': 'A',
        'G': 'C', 'C': 'G',
        'N': 'N'
    }
    
    sequence = sequence.upper()
    complement = ''.join(complement_map.get(base, 'N') for base in sequence)
    return complement[::-1]  # Reverse the string


def generate_random_dna(length: int, seed: int = None) -> str:
    """
    Generate a random DNA sequence.
    
    Args:
        length (int): Length of the sequence to generate
        seed (int): Random seed for reproducibility (optional)
    
    Returns:
        str: Random DNA sequence
    """
    import random
    if seed is not None:
        random.seed(seed)
    
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))


def extract_subsequence(sequence: str, start: int, end: int) -> str:
    """
    Extract a subsequence from a DNA sequence.
    
    Args:
        sequence (str): The full DNA sequence
        start (int): Starting position (0-indexed, inclusive)
        end (int): Ending position (0-indexed, exclusive)
    
    Returns:
        str: Extracted subsequence
    """
    return sequence[start:end]


def find_orfs(sequence: str, min_length: int = 100) -> List[Tuple[int, int, str]]:
    """
    Find Open Reading Frames (ORFs) in a DNA sequence.
    
    An ORF is a sequence that starts with ATG (start codon) and ends with
    a stop codon (TAA, TAG, or TGA).
    
    Args:
        sequence (str): DNA sequence
        min_length (int): Minimum ORF length in nucleotides (default: 100)
    
    Returns:
        List[Tuple[int, int, str]]: List of (start_pos, end_pos, orf_sequence)
    """
    sequence = sequence.upper()
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    
    # Search in all three reading frames
    for frame in range(3):
        i = frame
        while i < len(sequence) - 2:
            # Look for start codon
            if sequence[i:i+3] == start_codon:
                start_pos = i
                # Look for stop codon
                j = i + 3
                while j < len(sequence) - 2:
                    codon = sequence[j:j+3]
                    if codon in stop_codons:
                        end_pos = j + 3
                        orf_length = end_pos - start_pos
                        if orf_length >= min_length:
                            orfs.append((start_pos, end_pos, sequence[start_pos:end_pos]))
                        break
                    j += 3
                i = j  # Continue search after this ORF
            else:
                i += 3
    
    return orfs


def count_nucleotides(sequence: str) -> Dict[str, int]:
    """
    Count the occurrences of each nucleotide in a sequence.
    
    Args:
        sequence (str): DNA sequence
    
    Returns:
        Dict[str, int]: Dictionary mapping nucleotides to their counts
    """
    sequence = sequence.upper()
    return {
        'A': sequence.count('A'),
        'C': sequence.count('C'),
        'G': sequence.count('G'),
        'T': sequence.count('T'),
        'N': sequence.count('N')
    }


def split_sequence_into_chunks(sequence: str, chunk_size: int, overlap: int = 0) -> List[str]:
    """
    Split a long sequence into smaller chunks for processing.
    
    Args:
        sequence (str): DNA sequence to split
        chunk_size (int): Size of each chunk
        overlap (int): Number of bases to overlap between chunks (default: 0)
    
    Returns:
        List[str]: List of sequence chunks
    """
    chunks = []
    step = chunk_size - overlap
    
    for i in range(0, len(sequence), step):
        chunk = sequence[i:i + chunk_size]
        if len(chunk) >= chunk_size:  # Only include full chunks
            chunks.append(chunk)
        elif i == 0 or len(chunk) > overlap:  # Include last partial chunk if significant
            chunks.append(chunk)
    
    return chunks


def format_sequence_pretty(sequence: str, line_length: int = 80) -> str:
    """
    Format a DNA sequence for pretty printing.
    
    Args:
        sequence (str): DNA sequence
        line_length (int): Characters per line (default: 80)
    
    Returns:
        str: Formatted sequence with line breaks
    """
    lines = []
    for i in range(0, len(sequence), line_length):
        lines.append(sequence[i:i + line_length])
    return '\n'.join(lines)


def parse_fasta_contiguous(filepath: str) -> str:
    """Reads a FASTA/FNA file and returns the concatenated DNA sequence as a single string."""
    seq_parts = []
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip().upper())
    return "".join(seq_parts)
