#!/bin/bash

echo "=========================================="
echo "Generating 10 Synthetic DNA Datasets"
echo "=========================================="
echo ""

# Dataset 1: Small (100K bp)
echo "[1/10] Generating SYNTH_100K (100,000 bp)..."
python synthetic.py --length 100000 --seed 42 --name SYNTH_100K

# Dataset 2: Medium (250K bp)
echo ""
echo "[2/10] Generating SYNTH_250K (250,000 bp)..."
python synthetic.py --length 250000 --seed 43 --name SYNTH_250K

# Dataset 3: Medium-Large (500K bp)
echo ""
echo "[3/10] Generating SYNTH_500K (500,000 bp)..."
python synthetic.py --length 500000 --seed 44 --name SYNTH_500K

# Dataset 4: Large (750K bp)
echo ""
echo "[4/10] Generating SYNTH_750K (750,000 bp)..."
python synthetic.py --length 750000 --seed 45 --name SYNTH_750K

# Dataset 5: Very Large (1M bp)
echo ""
echo "[5/10] Generating SYNTH_1M (1,000,000 bp)..."
python synthetic.py --length 1000000 --seed 46 --name SYNTH_1M

# Dataset 6: Extra Large (2M bp)
echo ""
echo "[6/10] Generating SYNTH_2M (2,000,000 bp)..."
python synthetic.py --length 2000000 --seed 47 --name SYNTH_2M

# Dataset 7: GC-rich simulation (1M bp)
echo ""
echo "[7/10] Generating SYNTH_GC_RICH (1,000,000 bp)..."
python synthetic.py --length 1000000 --seed 48 --name SYNTH_GC_RICH

# Dataset 8: AT-rich simulation (1M bp)
echo ""
echo "[8/10] Generating SYNTH_AT_RICH (1,000,000 bp)..."
python synthetic.py --length 1000000 --seed 49 --name SYNTH_AT_RICH

# Dataset 9: Balanced (1M bp)
echo ""
echo "[9/10] Generating SYNTH_BALANCED (1,000,000 bp)..."
python synthetic.py --length 1000000 --seed 50 --name SYNTH_BALANCED

# Dataset 10: Repetitive patterns (1M bp)
echo ""
echo "[10/10] Generating SYNTH_REPETITIVE (1,000,000 bp)..."
python synthetic.py --length 1000000 --seed 51 --name SYNTH_REPETITIVE

echo ""
echo "=========================================="
echo "✓ All 10 datasets generated!"
echo "=========================================="
echo ""
echo "Verifying datasets..."
echo ""

# Count generated files
file_count=$(find synthetic_datasets -name "*.fna" | wc -l)
echo "Total .fna files: $file_count"

# List all datasets
echo ""
echo "Generated datasets:"
ls -lh synthetic_datasets/*/
