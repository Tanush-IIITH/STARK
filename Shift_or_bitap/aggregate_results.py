"""
Aggregate Benchmark Results for Shift-Or/Bitap Algorithm

This script scans all results folders and combines individual CSV files
into unified datasets for analysis.

Author: Shift-Or/Bitap Analysis
Date: November 2025
"""

import pandas as pd
import glob
import os
import json
from pathlib import Path

def aggregate_scaling_results():
    """
    Aggregate all scaling_results.csv files from:
    - results_exact/
    - results_approximate/
    - results_extended/
    - results_synthetic_exact/
    - results_synthetic_approximate/
    - results_synthetic_extended/
    """
    print("="*70)
    print("AGGREGATING SCALING RESULTS")
    print("="*70)

    all_data = []

    # Define result folders and their metadata
    folders = [
        ('results_exact', 'exact', None),
        ('results_approximate', 'approximate', None),
        ('results_extended', 'extended', None),
        ('results_synthetic_exact', 'exact', 'synthetic'),
        ('results_synthetic_approximate', 'approximate', 'synthetic'),
        ('results_synthetic_extended', 'extended', 'synthetic'),
    ]

    for folder, algo_type, data_source in folders:
        if not os.path.exists(folder):
            print(f"⚠ Skipping {folder} (not found)")
            continue

        # Find all scaling_results*.csv files
        if algo_type == 'approximate':
            # Handle k1, k2, k3 variants
            for k_val in [1, 2, 3]:
                pattern = f"{folder}/*/scaling_results_k{k_val}.csv"
                files = glob.glob(pattern)

                for filepath in files:
                    try:
                        df = pd.read_csv(filepath)
                        dataset_name = Path(filepath).parent.name

                        # Add metadata columns
                        df['dataset_name'] = dataset_name
                        df['algorithm_type'] = algo_type
                        df['k_value'] = k_val
                        df['data_source'] = data_source if data_source else 'real'
                        df['text_length_n'] = df['slice_size']  # Rename for consistency

                        all_data.append(df)
                    except Exception as e:
                        print(f"✗ Error reading {filepath}: {e}")
        else:
            # Exact and Extended
            pattern = f"{folder}/*/scaling_results.csv"
            files = glob.glob(pattern)

            for filepath in files:
                try:
                    df = pd.read_csv(filepath)
                    dataset_name = Path(filepath).parent.name

                    # Add metadata columns
                    df['dataset_name'] = dataset_name
                    df['algorithm_type'] = algo_type
                    df['k_value'] = None
                    df['data_source'] = data_source if data_source else 'real'
                    df['text_length_n'] = df['slice_size']

                    all_data.append(df)
                except Exception as e:
                    print(f"✗ Error reading {filepath}: {e}")

        print(f"✓ Processed {folder}: {len(glob.glob(f'{folder}/*/'))} datasets")

    # Combine all data
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)

        # Reorder columns
        cols = ['dataset_name', 'data_source', 'algorithm_type', 'k_value', 
                'text_length_n', 'avg_time_s', 'peak_memory_mb']
        available_cols = [c for c in cols if c in combined_df.columns]
        other_cols = [c for c in combined_df.columns if c not in available_cols]
        combined_df = combined_df[available_cols + other_cols]

        # Save
        output_file = 'benchmark_results_scaling.csv'
        combined_df.to_csv(output_file, index=False)
        print(f"\n✓ Saved {output_file}")
        print(f"  Total rows: {len(combined_df)}")
        print(f"  Datasets: {combined_df['dataset_name'].nunique()}")

        return combined_df
    else:
        print("\n✗ No scaling data found!")
        return None


def aggregate_pattern_length_results():
    """
    Aggregate all pattern_length_results.csv files.
    """
    print("\n" + "="*70)
    print("AGGREGATING PATTERN LENGTH RESULTS")
    print("="*70)

    all_data = []

    folders = [
        ('results_exact', 'exact', None),
        ('results_approximate', 'approximate', None),
        ('results_extended', 'extended', None),
        ('results_synthetic_exact', 'exact', 'synthetic'),
        ('results_synthetic_approximate', 'approximate', 'synthetic'),
        ('results_synthetic_extended', 'extended', 'synthetic'),
    ]

    for folder, algo_type, data_source in folders:
        if not os.path.exists(folder):
            continue

        if algo_type == 'approximate':
            # Handle k1, k2, k3 variants
            for k_val in [1, 2, 3]:
                pattern = f"{folder}/*/pattern_length_results_k{k_val}.csv"
                files = glob.glob(pattern)

                for filepath in files:
                    try:
                        df = pd.read_csv(filepath)
                        dataset_name = Path(filepath).parent.name

                        df['dataset_name'] = dataset_name
                        df['algorithm_type'] = algo_type
                        df['k_value'] = k_val
                        df['data_source'] = data_source if data_source else 'real'

                        all_data.append(df)
                    except Exception as e:
                        print(f"✗ Error reading {filepath}: {e}")
        else:
            pattern = f"{folder}/*/pattern_length_results.csv"
            files = glob.glob(pattern)

            for filepath in files:
                try:
                    df = pd.read_csv(filepath)
                    dataset_name = Path(filepath).parent.name

                    df['dataset_name'] = dataset_name
                    df['algorithm_type'] = algo_type
                    df['k_value'] = None
                    df['data_source'] = data_source if data_source else 'real'

                    all_data.append(df)
                except Exception as e:
                    print(f"✗ Error reading {filepath}: {e}")

        print(f"✓ Processed {folder}")

    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)

        # Reorder columns
        cols = ['dataset_name', 'data_source', 'algorithm_type', 'k_value', 
                'pattern_length', 'n', 'time_seconds_mean', 'peak_memory_mb']
        available_cols = [c for c in cols if c in combined_df.columns]
        other_cols = [c for c in combined_df.columns if c not in available_cols]
        combined_df = combined_df[available_cols + other_cols]

        output_file = 'benchmark_results_pattern.csv'
        combined_df.to_csv(output_file, index=False)
        print(f"\n✓ Saved {output_file}")
        print(f"  Total rows: {len(combined_df)}")
        print(f"  Datasets: {combined_df['dataset_name'].nunique()}")

        return combined_df
    else:
        print("\n✗ No pattern length data found!")
        return None


def generate_summary(scaling_df, pattern_df):
    """Generate summary statistics."""
    print("\n" + "="*70)
    print("GENERATING SUMMARY STATISTICS")
    print("="*70)

    summary = {
        'total_datasets': 0,
        'real_datasets': 0,
        'synthetic_datasets': 0,
        'algorithms': {},
        'generated_at': pd.Timestamp.now().isoformat()
    }

    if scaling_df is not None:
        summary['total_datasets'] = scaling_df['dataset_name'].nunique()
        summary['real_datasets'] = scaling_df[scaling_df['data_source'] == 'real']['dataset_name'].nunique()
        summary['synthetic_datasets'] = scaling_df[scaling_df['data_source'] == 'synthetic']['dataset_name'].nunique()

        for algo in scaling_df['algorithm_type'].unique():
            algo_data = scaling_df[scaling_df['algorithm_type'] == algo]
            summary['algorithms'][algo] = {
                'avg_time_s': float(algo_data['avg_time_s'].mean()),
                'avg_memory_mb': float(algo_data['peak_memory_mb'].mean()),
                'datasets_tested': int(algo_data['dataset_name'].nunique())
            }

    # Save summary
    with open('benchmark_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    print("\n✓ Saved benchmark_summary.json")
    print(f"\nSummary:")
    print(f"  Total datasets: {summary['total_datasets']}")
    print(f"  Real datasets: {summary['real_datasets']}")
    print(f"  Synthetic datasets: {summary['synthetic_datasets']}")

    return summary


def main():
    """Main execution."""
    print("\n" + "="*70)
    print("SHIFT-OR/BITAP: AGGREGATING BENCHMARK RESULTS")
    print("="*70 + "\n")

    # Aggregate scaling results
    scaling_df = aggregate_scaling_results()

    # Aggregate pattern length results
    pattern_df = aggregate_pattern_length_results()

    # Generate summary
    if scaling_df is not None or pattern_df is not None:
        summary = generate_summary(scaling_df, pattern_df)

    print("\n" + "="*70)
    print("✓ AGGREGATION COMPLETE!")
    print("="*70)
    print("\nGenerated files:")
    print("  - benchmark_results_scaling.csv")
    print("  - benchmark_results_pattern.csv")
    print("  - benchmark_summary.json")
    print("\nReady for analysis!")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
