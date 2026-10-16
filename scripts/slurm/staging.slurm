#!/bin/bash
# Run all staging tasks packed onto a single regular node
# using GNU Parallel instead of a SLURM array.
#
# Usage:
#   sbatch scripts/slurm/staging.slurm
#   PARAMS=scripts/params/staging_params_2700.tsv sbatch scripts/slurm/staging.slurm
#   sbatch --export=ALL,PARAMS=scripts/params/staging_params_2700.tsv scripts/slurm/staging.slurm
#
#SBATCH -J qpym2-staging
#SBATCH -A m2812
#SBATCH -q regular
#SBATCH -n 1
#SBATCH -C cpu
#SBATCH -t 04:00:00
#SBATCH -o /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_staging_%j.out
#SBATCH -e /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_staging_%j.err

set -euo pipefail

cd /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2

# Allow overriding the params TSV via env variable (e.g. for 2700 ROI run)
PARAMS="${PARAMS:-scripts/params/staging_params.tsv}"
JOBS_PARALLEL=32

N=$(awk -F'\t' 'NR>1 { c++ } END { print c+0 }' "${PARAMS}")
echo "Staging: ${N} tasks from ${PARAMS}, parallelism=${JOBS_PARALLEL}"

# Activate conda
source /global/common/software/cuore/sachi/qpym2/etc/profile.d/conda.sh
set +u
conda activate qpym2 || conda activate /global/common/software/cuore/sachi/qpym2
set -u

# Run all tasks from the TSV concurrently
awk -F'\t' 'NR>1 { print $NF }' "${PARAMS}" \
    | tr -d '\r' \
    | parallel -j "${JOBS_PARALLEL}" --halt soon,fail=10% \
        --joblog "logs/staging_parallel_${SLURM_JOB_ID}.log" \
        --line-buffer
