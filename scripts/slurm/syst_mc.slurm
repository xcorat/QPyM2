#!/bin/bash
# MC (BM-stat) systematic fit array job — one task per co_shift value.
# Fixed model point: co_cut=gamma10_v150, ndbd_cutoff=u1740.
#
# Generate params then submit:
#   python scripts/gen_params.py syst_mc
#   sbatch --array=1-7%7 scripts/slurm/syst_mc.slurm
#
# Or use the submit wrapper:
#   bash scripts/submit_syst.sh mc
#
#SBATCH -J qpym2-syst-mc
#SBATCH -A m2812
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH -t 04:00:00
#SBATCH --array=1-7%7
#SBATCH --image=docker:debian:bookworm-slim
#SBATCH -o /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_syst_mc_%A_%a.out
#SBATCH -e /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/logs/slurm_syst_mc_%A_%a.err

set -euo pipefail

export PARAMS="scripts/params/syst_mc_params.tsv"
cd /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2

srun shifter bash scripts/run_job.sh
