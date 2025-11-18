# Suffix Arrays & Trees Module

This directory houses the benchmarking code, datasets, and write-ups that back the suffix array/tree portion of the project. Use it as the single entry point for reproducing results or extending the experiments.

## Contents at a Glance
- `analysis.ipynb` – End-to-end notebook that loads the benchmark outputs, compares the algorithms across synthetic and NCBI datasets, and drives the bonus visualizations.
- `benchmark.py` – Runner that constructs each index, measures construction/query/memory metrics, and logs them to `benchmark_results.csv`.
- `benchmark_results.csv` – Cached measurements produced by `benchmark.py`; the notebook reads this file to plot the three primary graphs.
- `dataset/` – FASTA files (e.g., `GCA_…` genomes) that feed both the benchmarks and the interactive bonus cells.
- `graphs/` – All generated figures (`graph1_*`, `graph2_*`, `graph3_*`, `graph4_match_location_map.png`, suffix-tree PNGs, etc.).
- `suffix_array.py` – Pure-Python suffix array implementations (naive sorter plus Manber–Myers) along with helper routines such as pattern-location helpers.
- `suffix_tree.py` – Naive and Ukkonen suffix tree implementations plus traversal utilities that the Graphviz renderer relies on.
- `synthetic.py` – Utility script that emits synthetic DNA datasets in FASTA format to broaden the benchmark coverage.
- `SuffixTheoriticalAnalysis.pdf` – Companion document explaining the algorithms and their complexity proofs in a theoretical narrative.

Feel free to open the notebook first to explore the visualizations, then drop into the individual scripts if you need to recompute benchmarks or tweak the data structures.
