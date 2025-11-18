import os
import csv
import time
import psutil
import random
from shift_or_bitap import ShiftOr


# ------------------------------------------------------------
# MEMORY MEASUREMENT
# ------------------------------------------------------------
def current_memory_bytes():
    return psutil.Process(os.getpid()).memory_info().rss


# ------------------------------------------------------------
# FASTA READER
# ------------------------------------------------------------
def parse_fasta_contiguous(filepath: str) -> str:
    seq_parts = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip().upper())
    return "".join(seq_parts)


# ------------------------------------------------------------
# FIND FIRST DATASET FILE
# ------------------------------------------------------------
def find_first_dataset_file(dataset_root):
    patterns = ["**/*.fna", "**/*.fa", "**/*.fasta"]
    files = []
    for pat in patterns:
        files.extend(list(
            __import__("glob").glob(os.path.join(dataset_root, pat), recursive=True)
        ))

    files.sort(key=lambda x: os.path.getsize(x))
    return files[0] if files else None


# ------------------------------------------------------------
# SHIFT-OR CONSTRUCTION (with timing + memory)
# ------------------------------------------------------------
def so_construct_and_measure(pattern: str, max_errors: int = 0):
    start_mem = current_memory_bytes()
    start = time.perf_counter()

    try:
        so = ShiftOr(pattern, max_errors)
    except Exception:
        so = None

    elapsed = time.perf_counter() - start
    end_mem = current_memory_bytes()
    return so, elapsed, (end_mem - start_mem)


# ------------------------------------------------------------
# SHIFT-OR SEARCH TIME
# ------------------------------------------------------------
def so_search_time(so, text):
    start = time.perf_counter()
    matches = so.search(text)
    elapsed = time.perf_counter() - start
    return elapsed, matches


# ------------------------------------------------------------
# CSV HEADER
# ------------------------------------------------------------
def write_header(writer):
    writer.writerow([
        "algorithm",
        "dataset_name",
        "text_length_n",
        "pattern_length_m",
        "metric_type",
        "time_sec",
        "memory_bytes",
        "matches"
    ])


# ------------------------------------------------------------
# BENCHMARK 1: VARYING TEXT LENGTH
# ------------------------------------------------------------
def benchmark_varying_text_lengths(full_text, dataset_name, writer):

    fixed_pattern = "GAATTC"

    # text sizes (clamped to dataset length)
    n_values = [1_000, 2_000, 5_000, 10_000, 25_000, 50_000,
                75_000, 100_000, 150_000, 200_000]
    n_values = [n for n in n_values if n <= len(full_text)]

    # construct shift-or
    so, c_time, c_mem = so_construct_and_measure(fixed_pattern)

    writer.writerow(["shift_or", dataset_name, 0, len(fixed_pattern),
                     "construction", c_time, c_mem, 0])

    # search for each text length
    for n in n_values:
        text = full_text[:n]
        if so is not None:
            t_time, matches = so_search_time(so, text)
            writer.writerow(["shift_or", dataset_name, n, len(fixed_pattern),
                             "search", t_time, 0, len(matches)])
        else:
            # naive fallback
            start = time.perf_counter()
            matches = []
            idx = text.find(fixed_pattern)
            while idx != -1:
                matches.append(idx)
                idx = text.find(fixed_pattern, idx + 1)
            t_time = time.perf_counter() - start

            writer.writerow(["shift_or", dataset_name, n, len(fixed_pattern),
                             "search", t_time, 0, len(matches)])


# ------------------------------------------------------------
# BENCHMARK 2: VARYING PATTERN LENGTH
# ------------------------------------------------------------
def benchmark_varying_pattern_lengths(full_text, dataset_name, writer):

    n_fixed = min(100_000, len(full_text))
    text = full_text[:n_fixed]

    pattern_lengths = [5, 10, 20, 50, 100, 200]

    for m in pattern_lengths:

        if m >= n_fixed:
            continue

        # deterministic pattern location
        start_index = (n_fixed // 3) % (n_fixed - m)
        pattern = text[start_index:start_index + m]

        # construct shift-or
        so, c_time, c_mem = so_construct_and_measure(pattern)

        writer.writerow(["shift_or", dataset_name, n_fixed, m,
                         "construction", c_time, c_mem, 0])

        # search
        if so is not None:
            t_time, matches = so_search_time(so, text)
            writer.writerow(["shift_or", dataset_name, n_fixed, m,
                             "search", t_time, 0, len(matches)])

        else:
            # fallback naive search
            start = time.perf_counter()
            matches = []
            idx = text.find(pattern)
            while idx != -1:
                matches.append(idx)
                idx = text.find(pattern, idx + 1)
            t_time = time.perf_counter() - start

            writer.writerow(["shift_or", dataset_name, n_fixed, m,
                             "search", t_time, 0, len(matches)])


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
def main():

    # make dataset path relative to this directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_root = os.path.join(base_dir, "..", "DnA_dataset", "ncbi_dataset", "data")

    dataset_file = find_first_dataset_file(dataset_root)

    if not dataset_file:
        print(f"No dataset found under {dataset_root}.")
        return

    dataset_name = os.path.basename(dataset_file)
    print("Loading dataset:", dataset_name)

    full_text = parse_fasta_contiguous(dataset_file)

    print("Sequence length:", f"{len(full_text):,}", "bp")

    # write CSV
    with open("so_benchmark_results.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        write_header(writer)

        benchmark_varying_text_lengths(full_text, dataset_name, writer)
        benchmark_varying_pattern_lengths(full_text, dataset_name, writer)

    print("Benchmark completed. Output saved to so_benchmark_results.csv")


if __name__ == "__main__":
    main()
