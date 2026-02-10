#!/bin/bash
# scripts/run_parameter_sweeps.sh
#
# Run parameter sweep across all benchmark datasets to identify optimal
# clustering thresholds for haplotype discovery.
#
# Each sweep tests thresholds from 0.01 to 0.10 (or 0.08) in steps of 0.005.

set -e  # Exit on error

mkdir -p benchmarking/parameter_sweep

echo "Starting parameter sweeps across all datasets..."
echo "Each sweep will test 16-19 threshold values"
echo ""

# Threshold ranges (converted from min/max/step to comma-separated lists)
# 0.01 to 0.10 step 0.005 = 19 values
THRESHOLDS_10="0.010,0.015,0.020,0.025,0.030,0.035,0.040,0.045,0.050,0.055,0.060,0.065,0.070,0.075,0.080,0.085,0.090,0.095,0.100"
# 0.01 to 0.08 step 0.005 = 15 values
THRESHOLDS_8="0.010,0.015,0.020,0.025,0.030,0.035,0.040,0.045,0.050,0.055,0.060,0.065,0.070,0.075,0.080"

# Run sweeps (sequential by default; add --threads 4 to run phases in parallel)
echo "=== Sphyrnidae ==="
boldgenotyper-sweep data/Sphyrnidae.tsv \
    --thresholds $THRESHOLDS_10 \
    --output benchmarking/parameter_sweep/Sphyrnidae

echo ""
echo "=== Panulirus ==="
boldgenotyper-sweep data/Panulirus.tsv \
    --thresholds $THRESHOLDS_10 \
    --output benchmarking/parameter_sweep/Panulirus

echo ""
echo "=== Salmonidae ==="
boldgenotyper-sweep data/Salmonidae.tsv \
    --thresholds $THRESHOLDS_8 \
    --output benchmarking/parameter_sweep/Salmonidae

echo ""
echo "=== Pieridae ==="
boldgenotyper-sweep data/Pieridae.tsv \
    --thresholds $THRESHOLDS_10 \
    --output benchmarking/parameter_sweep/Pieridae

echo ""
echo "=== Carcharhiniformes ==="
boldgenotyper-sweep data/Carcharhiniformes.tsv \
    --thresholds $THRESHOLDS_8 \
    --output benchmarking/parameter_sweep/Carcharhiniformes

echo ""
echo "All parameter sweeps complete!"
echo "Results saved to: benchmarking/parameter_sweep/"
