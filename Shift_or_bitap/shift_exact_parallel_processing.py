#!/usr/bin/env python3
import os
import sys
import matplotlib
matplotlib.use('Agg')   # non-interactive backend (safe in multiprocessing)
from multiprocessing import Pool

# ---- import your existing benchmark function ----
from benchmark_exact import run_benchmark_on_dataset

# ---- Resolve paths relative to this script so no placeholders needed ----
script_dir = os.path.dirname(os.path.abspath(__file__))
DATASET_ROOT = os.path.abspath(os.path.join(script_dir, "../DnA_dataset/ncbi_dataset/data"))
# Put results next to your script (you can change this if you prefer another absolute path)
OUTPUT_ROOT = os.path.join(script_dir, "results_exact")

# Use all logical cores but leave 2 free: one for OS, one for your Jupyter notebook
NUM_CORES = max(1, os.cpu_count() - 2)

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
    os.makedirs(output_dir, exist_ok=True)

    # ---- IMPORTANT: switch cwd so relative writes go into output_dir ----
    prev_cwd = os.getcwd()
    try:
        os.chdir(output_dir)
        # call the existing function (unchanged)
        run_benchmark_on_dataset(fasta_file, output_dir)
        return f"Done: {dataset_name} -> {output_dir}"
    except Exception as e:
        return f"Error in {dataset_name}: {e}"
    finally:
        os.chdir(prev_cwd)


if __name__ == "__main__":
    # get sorted folder list
    all_datasets = sorted(
        d for d in os.listdir(DATASET_ROOT)
        if os.path.isdir(os.path.join(DATASET_ROOT, d))
    )

    print(f"Total datasets found: {len(all_datasets)} (root: {DATASET_ROOT})")

    # select folders 13..18 (1-indexed)
    datasets = all_datasets[12:18]
    print("Processing ONLY folders 13 to 18:")
    for name in datasets:
        print(" -", name)

    print(f"\nUsing {NUM_CORES} CPU cores...\n")
    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    with Pool(processes=NUM_CORES) as pool:
        results = pool.map(process_one_dataset, datasets)

    print("\n========== SUMMARY ==========")
    for r in results:
        print(r)
    print("=============================\n")
    print("Done.")


# import os
# import sys
# from multiprocessing import Pool
# from benchmark_exact import run_benchmark_on_dataset

# # ===========================================
# # CONFIGURATION
# # ===========================================
# DATASET_ROOT = "../DnA_dataset/ncbi_dataset/data"    # <-- UPDATE THIS
# OUTPUT_ROOT = "results_exact"

# # Use all CPU cores except 1 (to keep system responsive)
# NUM_CORES = max(1, os.cpu_count() - 1)

# # ===========================================
# # FUNCTION TO PROCESS ONE DATASET
# # ===========================================
# def process_one_dataset(dataset_name):
#     dataset_path = os.path.join(DATASET_ROOT, dataset_name)

#     # Find FASTA file inside dataset folder
#     fasta_file = None
#     for file in os.listdir(dataset_path):
#         if file.endswith(".fna") or file.endswith(".fasta"):
#             fasta_file = os.path.join(dataset_path, file)
#             break

#     if not fasta_file:
#         return f"Skipping {dataset_name} (no FASTA file)"

#     output_dir = os.path.join(OUTPUT_ROOT, dataset_name)

#     try:
#         run_benchmark_on_dataset(fasta_file, output_dir)
#         return f"Done: {dataset_name}"
#     except Exception as e:
#         return f"Error in {dataset_name}: {e}"


# # ===========================================
# # MAIN EXECUTION
# # ===========================================
# if __name__ == "__main__":

#     # Get all folders sorted
#     all_datasets = sorted([
#         d for d in os.listdir(DATASET_ROOT)
#         if os.path.isdir(os.path.join(DATASET_ROOT, d))
#     ])

#     print(f"Total datasets found: {len(all_datasets)}")

#     # Select only folders 13 → 18 (1-indexed)
#     # Python slicing: [12:18] gives index 12,13,14,15,16,17
#     datasets = all_datasets[13:17]

#     print(f"Running ONLY folders 14 to 17:")
#     for name in datasets:
#         print(" -", name)

#     print(f"\nUsing {NUM_CORES} CPU cores...\n")

#     os.makedirs(OUTPUT_ROOT, exist_ok=True)

#     # Run in parallel
#     with Pool(processes=NUM_CORES) as pool:
#         results = pool.map(process_one_dataset, datasets)

#     print("\n========== SUMMARY ==========")
#     for r in results:
#         print(r)
#     print("=============================\n")

#     print("Processing complete!")
