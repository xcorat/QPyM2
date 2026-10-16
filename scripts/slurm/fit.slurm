#!/bin/bash
# Bayesian fit job — supports both array and packed (GNU Parallel) modes.
#
# Array mode (submit.sh fit):
#   sbatch --array=1-63%20 scripts/slurm/fit.slurm
#
# Packed mode (submit_packed.sh fit):
#   sbatch scripts/slurm/fit.slurm
#   Runs ALL tasks on one node via GNU Parallel.
#
#SBATCH -J qpym2-fits
#SBATCH -A m2812
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --nodes=1
#SBATCH -t 01:00:00
#SBATCH --image=docker:debian:bookworm-slim
#SBATCH -o /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/data_2syst/logs/slurm_fit_%j.out
#SBATCH -e /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/data_2syst/logs/slurm_fit_%j.err

set -euo pipefail

export PARAMS="${PARAMS:-scripts/params/fit_params.tsv}"
cd /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # Array mode: single task
    srun shifter bash scripts/run_job.sh
else
    # Packed mode: all tasks via GNU Parallel.
    # 15 concurrent fits × ~4 PyMC worker processes each ≈ 60 cores used.
    MAX_PARALLEL=15

    source /global/common/software/cuore/sachi/qpym2/etc/profile.d/conda.sh
    set +u
    conda activate qpym2 || conda activate /global/common/software/cuore/sachi/qpym2
    set -u

    awk -F'\t' 'NR>1 { print $NF }' "${PARAMS}" \
        | tr -d '\r' \
        | parallel -j "${MAX_PARALLEL}" --line-buffer --joblog logs/fit_parallel.log \
            'export PYTENSOR_FLAGS="base_compiledir=${TMPDIR:-/tmp}/.pytensor_{%}"; eval {}'
fi
