#!/bin/bash
# MC (BM-stat) systematic fit array job — one task per co_shift value.
# Output goes to data_2syst/syst_results; logs to data_2syst/logs/.
#
# Generate params then submit:
#   python scripts/gen_params.py syst_mc --environment m2v5-1 --analysis tech_note \
#       --n 500 --output-dir data_2syst/syst_results \
#       --bias-csv data_2syst/toy_results/bias_merged.csv --xmin-bias 15 \
#       --save-traces --stat-only
#   sbatch --array=1-7%7 scripts/slurm/syst_mc_2syst.slurm
#
# Or use the submit wrapper:
#   bash scripts/submit_syst.sh mc --environment m2v5-1 --analysis tech_note \
#       --n 500 --output-dir data_2syst/syst_results \
#       --bias-csv data_2syst/toy_results/bias_merged.csv --xmin-bias 15 \
#       --save-traces --stat-only
#
#SBATCH -J qpym2-syst-mc-2syst
#SBATCH -A m2812
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH -t 04:00:00
#SBATCH --array=1-7%7
#SBATCH --image=docker:debian:bookworm-slim
#SBATCH -o /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/data_2syst/logs/slurm_syst_mc_%A_%a.out
#SBATCH -e /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2/data_2syst/logs/slurm_syst_mc_%A_%a.err

set -euo pipefail

export PARAMS="scripts/params/syst_mc_params.tsv"
cd /global/homes/x/xcorat/Analysis/M2/QPyM2_v5.1/QPyM2

srun shifter bash scripts/run_job.sh
