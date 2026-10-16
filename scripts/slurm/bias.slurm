#!/bin/bash
# Production bias fit array job.
#
# Array range and throttle are set by scripts/submit.sh bias (overrides --array below).
# To submit manually:
#   python scripts/gen_params.py bias
#   sbatch --array=1-100%50 scripts/slurm/bias.slurm
#
#SBATCH -J qpym2-bias
#SBATCH -A m2812
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH -t 01:00:00
#SBATCH --array=1-63
#SBATCH --image=docker:debian:bookworm-slim
#SBATCH -o /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/data_2syst/logs/slurm_bias_%A_%a.out
#SBATCH -e /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/data_2syst/logs/slurm_bias_%A_%a.err

set -euo pipefail

export PARAMS="${PARAMS:-scripts/params/bias_params.tsv}"
cd /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2

srun shifter bash scripts/run_job.sh
