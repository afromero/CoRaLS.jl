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

# The default project path matches the existing OSC layout. Override it with
# CORALS_PROJECT_DIR if the checkout is elsewhere.
set -euo pipefail

PROJECT_DIR="${CORALS_PROJECT_DIR:-${HOME}/BeattyLab/CoRaLS.jl}"
JULIA_BIN="${JULIA_BIN:-julia}"
TASK="${SLURM_ARRAY_TASK_ID:?This script must be submitted as a Slurm array.}"
VAR="${VAR:?Set VAR to ALT, ICE, or SLOPE.}"

[[ -f "$PROJECT_DIR/Project.toml" ]] || {
    echo "CoRaLS project not found at $PROJECT_DIR. Set CORALS_PROJECT_DIR." >&2
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
