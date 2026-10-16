#!/bin/bash
# Submit wrapper for systematic SLURM array jobs.
#
# Usage:
#   bash scripts/submit_syst.sh lss --lss-data-path <path>
#   bash scripts/submit_syst.sh lss --lss-data-path <path> --co-shift 0.8 1.0 1.2
#   bash scripts/submit_syst.sh lss --lss-data-path <path> --test
#   bash scripts/submit_syst.sh lss --lss-data-path <path> -- -q debug
#   bash scripts/submit_syst.sh mc
#   bash scripts/submit_syst.sh mc --n 500 --test
set -euo pipefail

LSS_DATA_PATH_DEFAULT='/global/homes/x/xcorat/Software/cuore-nersc-modern_v4.0.0/cuoremc/ares/data/Nature2021'

# ---- Parse systematic type (first positional arg) ------------------------
if [[ $# -lt 1 ]]; then
    echo "Usage: bash scripts/submit_syst.sh <syst_type> [gen_params flags] [-- sbatch flags]" >&2
    echo "  syst_type: lss mc" >&2
    exit 1
fi

SYST_TYPE="$1"; shift

case "${SYST_TYPE}" in
    lss) RUN_TYPE="syst_lss" ;;
    mc)  RUN_TYPE="syst_mc"  ;;
    *)
        echo "ERROR: unknown syst type '${SYST_TYPE}'. Must be one of: lss mc" >&2
        exit 1
        ;;
esac

# ---- Split remaining args: before '--' → gen_params, after → sbatch ------
GEN_ARGS=()
SBATCH_EXTRA=()
TEST_MODE=0
SEP=0
for arg in "$@"; do
    if [[ "$arg" == "--" ]]; then SEP=1; continue; fi
    if [[ $SEP -eq 0 ]]; then
        if [[ "$arg" == "--test" ]]; then
            TEST_MODE=1
        else
            GEN_ARGS+=("$arg")
        fi
    else
        SBATCH_EXTRA+=("$arg")
    fi
done

# Inject default --lss-data-path if not already supplied
if [[ "${RUN_TYPE}" == "syst_lss" ]]; then
    has_lss_path=0
    for arg in "${GEN_ARGS[@]}"; do
        [[ "$arg" == "--lss-data-path" ]] && has_lss_path=1
    done
    if [[ $has_lss_path -eq 0 ]]; then
        GEN_ARGS+=("--lss-data-path" "${LSS_DATA_PATH_DEFAULT}")
        echo "Using default --lss-data-path: ${LSS_DATA_PATH_DEFAULT}"
    fi
fi

cd "$(dirname "$0")/.."   # always run from project root

TSV="scripts/params/${RUN_TYPE}_params.tsv"
SLURM="scripts/slurm/${RUN_TYPE}.slurm"

CORES_PER_NODE="${CORES_PER_NODE:-256}"
CPUS_PER_TASK=$(grep -m1 '#SBATCH.*--cpus-per-task' "$SLURM" \
    | grep -oP '\d+') || CPUS_PER_TASK=256

# ---- Generate params -------------------------------------------------------
python scripts/gen_params.py "${RUN_TYPE}" "${GEN_ARGS[@]}"

N=$(( $(wc -l < "$TSV") - 1 ))
if [[ $N -lt 1 ]]; then
    echo "No tasks generated in ${TSV}" >&2
    exit 1
fi

if (( CPUS_PER_TASK >= CORES_PER_NODE )); then
    THROTTLE=${N}
else
    THROTTLE=$(( CORES_PER_NODE / CPUS_PER_TASK ))
fi

THROTTLE_VAR="${RUN_TYPE^^}_ARRAY_THROTTLE"
THROTTLE="${!THROTTLE_VAR:-${THROTTLE}}"

# ---- Test mode -------------------------------------------------------------
SBATCH_MODE_ARGS=()
ARRAY_SPEC="1-${N}%${THROTTLE}"
if [[ ${TEST_MODE} -eq 1 ]]; then
    NTEST=$(( N < 2 ? N : 2 ))
    ARRAY_SPEC="1-${NTEST}%${NTEST}"
    SBATCH_MODE_ARGS=(--qos=debug --time=00:30:00)
    echo "TEST mode: submitting ${NTEST} tasks to debug qos for 30 minutes"
fi

# ---- Submit ----------------------------------------------------------------
echo "Submitting ${RUN_TYPE}: --array=${ARRAY_SPEC}  (${SBATCH_EXTRA[*]:-no extra flags})"
sbatch --array="${ARRAY_SPEC}" "${SBATCH_EXTRA[@]}" "${SBATCH_MODE_ARGS[@]}" "$SLURM"
