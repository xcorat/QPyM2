#!/bin/bash
# Submit wrapper for SLURM *array* jobs (one task per array element).
# For packed run types that use GNU Parallel on a single node
# (input, staging), use submit_packed.sh instead.
#
# Usage:
#   bash scripts/submit.sh fit                                # main Bayesian fits
#   bash scripts/submit.sh toys                               # toy batch fits
#   bash scripts/submit.sh sens                               # sensitivity batch fits
#   bash scripts/submit.sh bias                               # bias batch fits
#   bash scripts/submit.sh toys --n-fits 100                  # smoke test
#   bash scripts/submit.sh bias --n-fits 50 -- -q debug      # override sbatch flags
#   bash scripts/submit.sh fit --no-regen                     # skip TSV regeneration, use existing params
set -euo pipefail

# ---- Parse run type (first positional arg) --------------------------------
if [[ $# -lt 1 ]]; then
    echo "Usage: bash scripts/submit.sh <run_type> [gen_params flags] [-- sbatch flags]" >&2
    echo "  run_type: fit | toys | sens | bias" >&2
    echo "  (for input/staging, use submit_packed.sh)" >&2
    exit 1
fi

RUN_TYPE="$1"; shift

# ---- Validate run type ----------------------------------------------------
case "${RUN_TYPE}" in
    input|staging)
        echo "ERROR: '${RUN_TYPE}' uses GNU Parallel (packed on one node), not SLURM arrays." >&2
        echo "  Use: bash scripts/submit_packed.sh ${RUN_TYPE} $*" >&2
        exit 1
        ;;
    fit|toys|sens|bias) ;;
    *)
        echo "ERROR: unknown run type '${RUN_TYPE}'. Must be one of: fit toys sens bias" >&2
        echo "  (for input/staging, use submit_packed.sh)" >&2
        exit 1
        ;;
esac

# ---- Split remaining args: before '--' → gen_params, after → sbatch ------
GEN_ARGS=()
SBATCH_EXTRA=()
TEST_MODE=0
NO_REGEN=0
SEP=0
for arg in "$@"; do
    if [[ "$arg" == "--" ]]; then SEP=1; continue; fi
    if [[ $SEP -eq 0 ]]; then
        if [[ "$arg" == "--test" ]]; then
            TEST_MODE=1
        elif [[ "$arg" == "--no-regen" ]]; then
            NO_REGEN=1
        else
            GEN_ARGS+=("$arg")
        fi
    else
        SBATCH_EXTRA+=("$arg")
    fi
done

cd "$(dirname "$0")/.."   # always run from project root

# ---- Per-type defaults ----------------------------------------------------
TSV="scripts/params/${RUN_TYPE}_params.tsv"
SLURM="scripts/slurm/${RUN_TYPE}.slurm"

# Perlmutter CPU: 128 physical cores, 256 hardware threads (SMT-2).
CORES_PER_NODE="${CORES_PER_NODE:-256}"

# Read --cpus-per-task from the SLURM spec (fallback: 2).
CPUS_PER_TASK=$(grep -m1 '#SBATCH.*--cpus-per-task' "$SLURM" \
    | grep -oP '\d+') || CPUS_PER_TASK=2

# ---- Generate params ------------------------------------------------------
if [[ ${NO_REGEN} -eq 1 ]]; then
    echo "--no-regen: skipping gen_params, using existing ${TSV}"
else
    python scripts/gen_params.py "${RUN_TYPE}" "${GEN_ARGS[@]}"
fi

N=$(( $(wc -l < "$TSV") - 1 ))   # exclude header
if [[ $N -lt 1 ]]; then
    echo "No tasks generated in ${TSV}" >&2
    exit 1
fi

# When each task requests a full node, run all tasks concurrently (throttle=N).
if (( CPUS_PER_TASK >= CORES_PER_NODE )); then
    THROTTLE=${N}
else
    THROTTLE=$(( CORES_PER_NODE / CPUS_PER_TASK ))
fi

# Allow per-type env-var override (e.g. TOYS_ARRAY_THROTTLE=10).
THROTTLE_VAR="${RUN_TYPE^^}_ARRAY_THROTTLE"
THROTTLE="${!THROTTLE_VAR:-${THROTTLE}}"

# ---- Test mode (staging-style debug queue) ---------------------------------
SBATCH_MODE_ARGS=()
ARRAY_SPEC="1-${N}%${THROTTLE}"
if [[ ${TEST_MODE} -eq 1 ]]; then
    NTEST=$(( N < 2 ? N : 2 ))
    ARRAY_SPEC="1-${NTEST}%${NTEST}"
    SBATCH_MODE_ARGS=(--qos=debug --time=00:10:00)
    echo "TEST mode: submitting ${NTEST} tasks to debug qos for 10 minutes"
fi

# ---- Submit ---------------------------------------------------------------
echo "Submitting --array=${ARRAY_SPEC}  (${SBATCH_EXTRA[*]:-no extra flags})"
sbatch --array="${ARRAY_SPEC}" "${SBATCH_EXTRA[@]}" "${SBATCH_MODE_ARGS[@]}" "$SLURM"
