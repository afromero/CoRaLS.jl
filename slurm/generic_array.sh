#!/usr/bin/env bash
#SBATCH -A PAS0654
#SBATCH --job-name=accpt_array
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=slurm/out/corals_%A_%a.out
#SBATCH --error=slurm/err/corals_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=100:30:00
##SBATCH --partition=hugemem
#SBATCH --mem=50G

# The default project path matches Machtay's OSC layout. Override
# CORALS_PROJECT_DIR to select a sibling worktree for a particular job.
set -euo pipefail

CORALS_OSC_ROOT="${CORALS_OSC_ROOT:-/fs/project/PAS0654/machtay_corals}"
PROJECT_DIR="${CORALS_PROJECT_DIR:-${CORALS_OSC_ROOT}/CoRaLS.jl}"
CORALS_JULIA_DEPOT="${CORALS_JULIA_DEPOT:-${CORALS_OSC_ROOT}/julia_depot}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-${CORALS_JULIA_DEPOT}:}"
TASK="${SLURM_ARRAY_TASK_ID:?This script must be submitted as a Slurm array.}"
VAR="${VAR:?Set VAR to ALT, ICE, or SLOPE.}"

# Batch shells are often non-interactive, so they do not necessarily inherit
# the Julia module loaded in the submitting shell.  An explicit JULIA_BIN
# remains authoritative; otherwise load OSC's Julia 1.11 module here.
if [[ -z "${JULIA_BIN:-}" ]] && ! command -v julia >/dev/null 2>&1; then
    if ! type module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
        # shellcheck disable=SC1091
        source /etc/profile.d/modules.sh
    fi
    type module >/dev/null 2>&1 && module load julia/1.11
fi
JULIA_BIN="${JULIA_BIN:-julia}"

[[ -f "$PROJECT_DIR/Project.toml" ]] || {
    echo "CoRaLS project not found at $PROJECT_DIR. Set CORALS_PROJECT_DIR." >&2
    exit 2
}
command -v "$JULIA_BIN" >/dev/null 2>&1 || {
    echo "Julia executable not found: $JULIA_BIN" >&2
    echo "Load julia/1.11, or submit with JULIA_BIN=/full/path/to/julia." >&2
    exit 2
}

case "$VAR" in
    ALT)
        ALT="$TASK"
        ;;
    ICE)
        ICE="$TASK"
        ;;
    SLOPE)
        : "${SLOPE_MODELS:?Set SLOPE_MODELS to a colon-separated ordered list.}"
        [[ "$TASK" =~ ^[0-9]+$ ]] && (( TASK >= 1 )) || {
            echo "SLOPE array task ID must be a positive integer; got '$TASK'." >&2
            exit 2
        }
        IFS=':' read -r -a slope_models <<< "$SLOPE_MODELS"
        (( TASK <= ${#slope_models[@]} )) || {
            echo "Task $TASK has no model in SLOPE_MODELS='$SLOPE_MODELS'." >&2
            exit 2
        }
        SLOPE_MODEL="${slope_models[TASK - 1]}"
        ;;
    *)
        echo "VAR must be ALT, ICE, or SLOPE; got '$VAR'." >&2
        exit 2
        ;;
esac

: "${ALT:?Set ALT to the altitude index (the Julia script multiplies it by 5 km).}"
: "${ENERGY:?Set ENERGY to the energy multiplier.}"
: "${ICE:?Set ICE to the ice depth in metres.}"
: "${ANT:?Set ANT to the number of antennas.}"
: "${TRIG:?Set TRIG to the trigger multiplicity.}"
: "${ANG:?Set ANG to the antenna pointing angle in degrees.}"
: "${FREQ1:?Set FREQ1 to the minimum frequency in MHz.}"
: "${TEXP:?Set TEXP to log10(trials per bin before the altitude multiplier).}"
SLOPE_MODEL="${SLOPE_MODEL:-gaussian_0}"

mkdir -p slurm/out slurm/err
cd "$PROJECT_DIR"

echo "alt=$((5 * ALT)) km  energyMult=$ENERGY  ice=$ICE m  ant=$ANT  trig=$TRIG  angle=$ANG deg  freqMin=$FREQ1 MHz  TEXP=$TEXP  slope=$SLOPE_MODEL"

"$JULIA_BIN" --project="$PROJECT_DIR" "$PROJECT_DIR/slurm/generic_acceptance.jl" \
    "$ALT" "$ENERGY" "$ICE" "$ANT" "$TRIG" "$ANG" "$FREQ1" "$TEXP" "$SLOPE_MODEL"
