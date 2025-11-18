import os
import sys
import matplotlib
matplotlib.use('Agg')

from multiprocessing import Pool
from benchmark_exact import run_benchmark_on_dataset

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATASET_ROOT = os.path.join(BASE_DIR, "synthetic_datasets")
OUTPUT_ROOT = os.path.join(BASE_DIR, "results_synthetic_exact")

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

    # Output folder
    output_dir = os.path.join(OUTPUT_ROOT, dataset_name)
    os.makedirs(output_dir, exist_ok=True)

    # Run benchmark
    old = os.getcwd()
    os.chdir(output_dir)

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

    print("Running EXACT for synthetic datasets:")
    for d in all_folders:
        print(" -", d)

    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    with Pool(NUM_CORES) as p:
        results = p.map(process_one_dataset, all_folders)

    print("\n===== SUMMARY =====")
    for r in results:
        print(r)
