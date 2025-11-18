import os
import sys
from multiprocessing import Pool
from benchmark_extended import run_benchmark_on_dataset

# ===============================
# CONFIGURATION
# ===============================

DATASET_ROOT = "../DnA_dataset/ncbi_dataset/data"   # <-- UPDATE THIS
OUTPUT_ROOT = "results_extended"

# Number of CPU cores to use
NUM_CORES = max(1, os.cpu_count() - 2)   # uses (cores - 1)

# ===============================
# FUNCTION FOR ONE DATASET
# ===============================
def process_one_dataset(dataset_name):
    dataset_path = os.path.join(DATASET_ROOT, dataset_name)

    # Find FASTA file inside dataset folder
    fasta_file = None
    for file in os.listdir(dataset_path):
        if file.endswith(".fna") or file.endswith(".fasta"):
            fasta_file = os.path.join(dataset_path, file)
            break

    if not fasta_file:
        return f"Skipping {dataset_name} (no FASTA file)"

    output_dir = os.path.join(OUTPUT_ROOT, dataset_name)

    try:
        run_benchmark_on_dataset(fasta_file, output_dir)
        return f"Done: {dataset_name}"
    except Exception as e:
        return f"Error in {dataset_name}: {e}"


# ===============================
# MAIN PARALLEL EXECUTION
# ===============================
if __name__ == "__main__":

    # find all dataset folders
    datasets = sorted([
        d for d in os.listdir(DATASET_ROOT)
        if os.path.isdir(os.path.join(DATASET_ROOT, d))
    ])

    print(f"Found {len(datasets)} datasets")
    print(f"Running in parallel using {NUM_CORES} CPU cores...\n")

    # Create results directory if missing
    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    # Run all datasets in parallel
    with Pool(processes=NUM_CORES) as pool:
        results = pool.map(process_one_dataset, datasets)

    print("\n========= SUMMARY =========")
    for r in results:
        print(r)
    print("===========================")

    print("\nAll tasks completed!")
