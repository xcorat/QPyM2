#!/bin/bash
# Generic job runner — called by shifter inside each SLURM array task.
# Looks up the command for this task in a TSV file by SLURM_ARRAY_TASK_ID,
# activates the conda env, then executes the command.
#
# The PARAMS env var must be set by the calling SLURM script
# (scripts/slurm/*.slurm each set it explicitly).
# Fallback default is toys_params.tsv.
set -euo pipefail

PARAMS="${PARAMS:-scripts/params/dummy.tsv}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set" >&2
    exit 1
fi

CMD=$(awk -F'\t' -v id="${SLURM_ARRAY_TASK_ID}" \
    'NR>1 && $1==id { print $NF; exit }' "${PARAMS}")
CMD="${CMD//$'\r'/}"

if [[ -z "${CMD}" ]]; then
    echo "ERROR: no row found for task_id=${SLURM_ARRAY_TASK_ID} in ${PARAMS}" >&2
    exit 1
fi

echo "[task ${SLURM_ARRAY_TASK_ID}] ${CMD}"

source /global/common/software/cuore/sachi/qpym2/etc/profile.d/conda.sh
set +u
conda activate qpym2 || conda activate /global/common/software/cuore/sachi/qpym2
set -u

PYTENSOR_BASE="${TMPDIR:-/tmp}/.pytensor_${SLURM_JOB_ID:-nojid}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${PYTENSOR_BASE}"
export PYTENSOR_FLAGS="base_compiledir=${PYTENSOR_BASE}"

echo "Running command in $(pwd): ${CMD}"
T_START=$SECONDS
eval "${CMD}"
echo "[task ${SLURM_ARRAY_TASK_ID}] elapsed: $((SECONDS - T_START))s"
