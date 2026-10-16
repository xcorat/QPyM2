#!/bin/bash
# Production toy fit array job.
#
# Array range and throttle are set by scripts/submit.sh toys (overrides --array below).
# To submit manually:
#   python scripts/gen_params.py toys
#   sbatch --array=1-14%14 scripts/slurm/toys.slurm
#
#SBATCH -J qpym2-toys
#SBATCH -A m2812
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH -t 01:00:00
#SBATCH --array=1-14
#SBATCH --image=docker:debian:bookworm-slim
#SBATCH -o /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_toys_%A_%a.out
#SBATCH -e /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_toys_%A_%a.err

set -euo pipefail

export PARAMS="${PARAMS:-scripts/params/toys_params.tsv}"
cd /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2

srun shifter bash scripts/run_job.sh
