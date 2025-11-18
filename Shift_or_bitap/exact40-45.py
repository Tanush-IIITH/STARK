import os
import sys
import matplotlib
matplotlib.use('Agg')

from multiprocessing import Pool
from benchmark_exact import run_benchmark_on_dataset

DATASET_ROOT = "/home/keshav-goel/Desktop/AAD/STARK/DnA_dataset/ncbi_dataset/data"
OUTPUT_ROOT = "/home/keshav-goel/Desktop/AAD/STARK/Shift_or_bitap/results_exact"
NUM_CORES = max(1, os.cpu_count() - 1)


def process_one_dataset(dataset_name):

    dataset_path = os.path.join(DATASET_ROOT, dataset_name)

    # Find FASTA
    fasta_file = None
    for f in os.listdir(dataset_path):
        if f.endswith(".fna") or f.endswith(".fasta"):
            fasta_file = os.path.join(dataset_path, f)
            break

    if not fasta_file:
        return f"Skipping {dataset_name} (no FASTA)"

    output_dir = os.path.join(OUTPUT_ROOT, dataset_name)
    os.makedirs(output_dir, exist_ok=True)

    old = os.getcwd()
    os.chdir(output_dir)
    
    from shift_or_utils import read_fasta_file
    seqs = read_fasta_file(fasta_file)

    if seqs:
        key = list(seqs.keys())[0]
        text = seqs[key]

        LIMIT = 100000   # or 10000 or 1000 — YOU decide
        if len(text) > LIMIT:
            print(f"  ⚠ Truncating genome from {len(text)} bp to {LIMIT} bp")
            # overwrite sequence with truncated version
            seqs[key] = text[:LIMIT]

        # Monkey patch the benchmark’s read_fasta_file to return this truncated seq
        import benchmark_exact
        benchmark_exact.read_fasta_file = lambda _: seqs


    try:
        run_benchmark_on_dataset(fasta_file, output_dir)
        return f"Done: {dataset_name}"
    except Exception as e:
        return f"Error in {dataset_name}: {e}"
    finally:
        os.chdir(old)


if __name__ == "__main__":
    all_folders = sorted(
        d for d in os.listdir(DATASET_ROOT)
        if os.path.isdir(os.path.join(DATASET_ROOT, d))
    )

    datasets = all_folders[39:45]  # folders 40 → 45

    print("Running EXACT for datasets 40 → 45:")
    for d in datasets:
        print(" -", d)

    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    with Pool(NUM_CORES) as p:
        results = p.map(process_one_dataset, datasets)

    print("\n===== SUMMARY =====")
    for r in results:
        print(r)
