#!/bin/bash
# Run all input-creation tasks packed onto a single regular node
# using GNU Parallel instead of a SLURM array.
#
# Each task handles one (co_cut, ndbd_cutoff) pair and loops over all
# co_shifts internally, caching shared ROOT reads across shifts.
# Running all tasks on the same node lets the OS page cache serve
# staging ROOT files from RAM after the first process reads them.
#
# Usage:
#   python scripts/gen_params.py input
#   sbatch scripts/slurm/input.slurm
#
#SBATCH -J qpym2-input
#SBATCH -A m2812
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --nodes=1
#SBATCH -t 05:45:00
#SBATCH -o /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_input_%j.out
#SBATCH -e /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_input_%j.err

set -euo pipefail

cd /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2

PARAMS="scripts/params/input_params.tsv"
N_JOBS=$(awk -F'\t' 'NR>1 { c++ } END { print c }' "${PARAMS}")

# Activate conda
source /global/common/software/cuore/sachi/qpym2/etc/profile.d/conda.sh
set +u
conda activate qpym2 || conda activate /global/common/software/cuore/sachi/qpym2
set -u

# Run all tasks from the TSV concurrently on this node
awk -F'\t' 'NR>1 { print $NF }' "${PARAMS}" \
    | tr -d '\r' \
    | parallel -j "${N_JOBS}" --line-buffer --halt soon,fail=10% --joblog logs/input_parallel.log
