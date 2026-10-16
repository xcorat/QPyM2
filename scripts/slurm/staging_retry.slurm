#!/bin/bash
# Re-run failed staging tasks (33-86) packed onto a single node
# using GNU Parallel instead of a SLURM array.
#
# Usage:
#   sbatch scripts/slurm/staging_retry.slurm
#
#SBATCH -J qpym2-staging-retry
#SBATCH -A m2812
#SBATCH -q regular
#SBATCH -n 1
#SBATCH -C cpu
#SBATCH -t 03:10:00
#SBATCH --image=docker:debian:bookworm-slim
#SBATCH -o /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_staging_retry_%j.out
#SBATCH -e /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_staging_retry_%j.err

set -euo pipefail

cd /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2

PARAMS="scripts/params/staging_params.tsv"
TASK_START=1
TASK_END=86
# 54 failed tasks (33..86), all run concurrently on the node
JOBS_PARALLEL=32

# Activate conda
source /global/common/software/cuore/sachi/qpym2/etc/profile.d/conda.sh
set +u
conda activate qpym2 || conda activate /global/common/software/cuore/sachi/qpym2
set -u

# Extract commands for failed tasks (task_id 33..86) from the TSV
awk -F'\t' -v lo="${TASK_START}" -v hi="${TASK_END}" \
    'NR>1 && $1>=lo && $1<=hi { print $NF }' "${PARAMS}" \
    | parallel -j "${JOBS_PARALLEL}" --halt soon,fail=10% --joblog logs/staging_retry_parallel.log
