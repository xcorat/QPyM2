#!/bin/bash
# Submit wrapper for "packed" run types that use GNU Parallel on a single node.
# These SLURM scripts run ALL tasks from the TSV internally via `parallel`,
# so they must NOT be submitted as SLURM array jobs.
#
# Supported run types: input, staging
#
# Usage:
#   bash scripts/submit_packed.sh input                            # input HDF5 creation
#   bash scripts/submit_packed.sh input --analysis tech_note       # with gen_params flags
#   bash scripts/submit_packed.sh input --test                     # debug queue smoke test
#   bash scripts/submit_packed.sh staging                          # staging ROOT creation
#   bash scripts/submit_packed.sh staging -- -t 04:00:00           # override sbatch flags
set -euo pipefail

# ---- Parse run type (first positional arg) --------------------------------
if [[ $# -lt 1 ]]; then
    echo "Usage: bash scripts/submit_packed.sh <run_type> [gen_params flags] [-- sbatch flags]" >&2
    echo "  run_type: input | staging | fit | toys" >&2
    exit 1
fi

RUN_TYPE="$1"; shift

# ---- Validate run type ----------------------------------------------------
case "${RUN_TYPE}" in
    input|staging|fit|toys) ;;
    *)
        echo "ERROR: '${RUN_TYPE}' is not a valid run type." >&2
        echo "  Supported: input, staging, fit, toys" >&2
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

if [[ ! -f "$SLURM" ]]; then
    echo "ERROR: SLURM script not found: ${SLURM}" >&2
    exit 1
fi

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

# ---- Test mode (debug queue) ----------------------------------------------
SBATCH_MODE_ARGS=()
if [[ ${TEST_MODE} -eq 1 ]]; then
    SBATCH_MODE_ARGS=(--qos=debug --time=00:10:00)
    echo "TEST mode: submitting to debug qos for 10 minutes"
fi

# ---- Submit (NO --array — the slurm script runs all tasks via parallel) ---
echo "Submitting ${RUN_TYPE}: ${N} tasks packed on 1 node  (${SBATCH_EXTRA[*]:-no extra flags})"
sbatch "${SBATCH_EXTRA[@]}" "${SBATCH_MODE_ARGS[@]}" "$SLURM"
