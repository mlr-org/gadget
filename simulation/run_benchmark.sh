#!/bin/bash
# Run full GADGET vs effector efficiency benchmark.
# Usage: bash simulation/run_benchmark.sh [small|large]
#
# Presets
#   small (default): N=500,1000,5000  D=5,10,20  fixed_N=1000  fixed_D=10
#   large:           N=5000,10000,25000,50000  D=10,20,50,100  fixed_N=10000  fixed_D=20
#
# Output locations (large uses separate dirs so a long run does not overwrite small results):
#   small: simulation/data/benchmark/, simulation/results/benchmark/, simulation/results/figures/
#   large: simulation/data/benchmark_large/, simulation/results/benchmark_large/, simulation/results/figures_large/
#
# Dependencies
#   R:      ranger, gadget (devtools::load_all), data.table, ggplot2, gridExtra, scales
#   Python: pip install effector numpy pandas scikit-learn
#   Optional: PYTHON=/path/to/python bash simulation/run_benchmark.sh large

set -e
cd "$(dirname "$0")/.."
PYTHON="${PYTHON:-python3}"

PRESET="${1:-small}"

if [ "$PRESET" = "large" ]; then
  N_VEC="5000,10000,25000,50000"
  D_VEC="10,20,50,100"
  FIXED_N=10000
  FIXED_D=20
  DATADIR="simulation/data/benchmark_large"
  OUTDIR="simulation/results/benchmark_large"
  FIGDIR="simulation/results/figures_large"
else
  N_VEC="500,1000,5000"
  D_VEC="5,10,20"
  FIXED_N=1000
  FIXED_D=10
  DATADIR="simulation/data/benchmark"
  OUTDIR="simulation/results/benchmark"
  FIGDIR="simulation/results/figures"
fi

N_GRID_VEC="10,20,50"
N_INT_VEC="10,20,50"
REPS=5
PREDICT_REPS=20

echo "=== Benchmark preset: ${PRESET} ==="
echo "    N_VEC=${N_VEC}  D_VEC=${D_VEC}  fixed_N=${FIXED_N}  fixed_D=${FIXED_D}"
echo "    DATADIR=${DATADIR}  OUTDIR=${OUTDIR}  FIGDIR=${FIGDIR}"

echo "1. Generating benchmark data..."
Rscript simulation/generate_benchmark_data.R \
  --outdir "${DATADIR}" \
  --N-vec "${N_VEC}" --D-vec "${D_VEC}"

echo "2. Running GADGET benchmark (R)..."
Rscript simulation/benchmark_gadget.R \
  --datadir "${DATADIR}" --outdir "${OUTDIR}" \
  --reps "${REPS}" --predict-reps "${PREDICT_REPS}" \
  --N-vec "${N_VEC}" --D-vec "${D_VEC}" \
  --n-grid-vec "${N_GRID_VEC}" --n-int-vec "${N_INT_VEC}" \
  --fixed-N "${FIXED_N}" --fixed-D "${FIXED_D}"

echo "3. Running effector benchmark (Python)..."
"${PYTHON}" -u simulation/benchmark_effector.py \
  --datadir "${DATADIR}" --outdir "${OUTDIR}" \
  --reps "${REPS}" --predict-reps "${PREDICT_REPS}" \
  --N-vec "${N_VEC}" --D-vec "${D_VEC}" \
  --n-grid-vec "${N_GRID_VEC}" --n-int-vec "${N_INT_VEC}" \
  --fixed-N "${FIXED_N}" --fixed-D "${FIXED_D}"

echo "4. Summarizing and plotting..."
Rscript simulation/summarize_benchmark.R \
  --indir "${OUTDIR}" --figdir "${FIGDIR}" \
  --fixed-N "${FIXED_N}" --fixed-D "${FIXED_D}"

echo "Done. Raw CSVs and summary.csv in ${OUTDIR}/; figures in ${FIGDIR}/"
