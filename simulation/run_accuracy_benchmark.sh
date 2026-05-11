#!/bin/bash
# Run full accuracy pipeline (plan §9): shared data, GADGET, effector, summarize.
# Usage: bash simulation/run_accuracy_benchmark.sh [--n-seeds 30] [--variants num_0,num_05,cat] ...
# Example (fast test):
#   bash simulation/run_accuracy_benchmark.sh --n-seeds 2 --N-vec 500 --D-vec 10 --variants num_0
#
# Requires: R (gadget load + data.table + ggplot2); Python effector, pandas, numpy.

set -e
cd "$(dirname "$0")/.."

echo "1. Generating shared accuracy datasets..."
Rscript simulation/generate_accuracy_data.R "$@"

echo "2. GADGET accuracy..."
Rscript simulation/benchmark_accuracy_gadget.R "$@"

echo "3. effector accuracy..."
python3 -u simulation/benchmark_accuracy_effector.py "$@"

echo "4. Summarize..."
Rscript simulation/summarize_accuracy.R

echo "Done. Results in simulation/results/accuracy/ and simulation/results/figures/accuracy_*.png"
