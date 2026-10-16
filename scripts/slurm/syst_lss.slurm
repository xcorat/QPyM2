#!/bin/bash
# LSS systematic fit array job (one task per co_shift).
#
# Each task runs all datasets for a single co_shift value using
# the runner's internal multiprocessing pool.
#
# Array range and throttle are set by scripts/submit_syst.sh lss.
# To submit manually:
#   python scripts/gen_params.py syst_lss --syst-type sscl
#   sbatch --array=1-7 scripts/slurm/syst_lss.slurm
#
#SBATCH -J qpym2-lss
#SBATCH -A m2812
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH -t 04:00:00
#SBATCH --array=1-7
#SBATCH --image=docker:debian:bookworm-slim
#SBATCH -o /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_lss_%A_%a.out
#SBATCH -e /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_lss_%A_%a.err

set -euo pipefail

export PARAMS="scripts/params/syst_lss_params.tsv"
cd /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2

srun shifter bash scripts/run_job.sh
