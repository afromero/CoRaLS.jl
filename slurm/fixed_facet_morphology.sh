#!/usr/bin/env bash
#SBATCH -A PAS0654
#SBATCH --job-name=facet_morphology
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm/out/facet_morphology_%j.out
#SBATCH --error=slurm/err/facet_morphology_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=16G

# Run a paired fixed-facet morphology scan on one OSC compute node.
#
# This job intentionally remains single-process: every tilt must evaluate the
# same candidate event, and the current Julia diagnostic is serial.  The
# recorder streams event batches to CSV, so the 16 GB request is for Julia,
# package compilation, and modest buffers rather than the full output table.
#
# Submit from the project root, for example:
#   mkdir -p slurm/out slurm/err
#   sbatch --export=ALL,CORALS_PROJECT_DIR="$PWD",MORPH_NTRIALS=1000000 \
#     slurm/fixed_facet_morphology.sh
#
# Optional exported settings:
#   MORPH_NTRIALS       shared candidates (default 500000)
#   MORPH_OUTPUT_DIR    relative or absolute output directory
#   MORPH_TILTS         comma-separated degrees (default 0,...,10)
#   MORPH_THETA0S       comma-separated LPDA boresights (default -48)
#   MORPH_WRITE_BATCH   candidates per CSV flush (default 5000)
#   MORPH_SKIP_PLOTS    true/1 to skip matplotlib figures after summaries
#   JULIA_BIN           Julia executable (default julia from module/PATH)

set -euo pipefail

CORALS_OSC_ROOT="${CORALS_OSC_ROOT:-/fs/project/PAS0654/machtay_corals}"
PROJECT_DIR="${CORALS_PROJECT_DIR:-${CORALS_OSC_ROOT}/CoRaLS.jl}"
CORALS_JULIA_DEPOT="${CORALS_JULIA_DEPOT:-${CORALS_OSC_ROOT}/julia_depot}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-${CORALS_JULIA_DEPOT}:}"
MORPH_NTRIALS="${MORPH_NTRIALS:-500000}"
MORPH_OUTPUT_DIR="${MORPH_OUTPUT_DIR:-data/terrain_diagnostics/fixed_facet_0to10_1deg_osc}"
MORPH_TILTS="${MORPH_TILTS:-0,1,2,3,4,5,6,7,8,9,10}"
MORPH_THETA0S="${MORPH_THETA0S:--48}"
MORPH_WRITE_BATCH="${MORPH_WRITE_BATCH:-5000}"
MORPH_SKIP_PLOTS="${MORPH_SKIP_PLOTS:-false}"

# Ensure compute-node batch shells have the tools needed by the Julia recorder
# and its Python cache/plot stages.  Explicit JULIA_BIN still takes precedence.
if ! type module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi
if type module >/dev/null 2>&1; then
    [[ -n "${JULIA_BIN:-}" ]] || module load julia/1.11
    module load python/3.12
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
command -v python3 >/dev/null 2>&1 || {
    echo "python3 is not on PATH after attempting to load python/3.12." >&2
    exit 2
}

cd "$PROJECT_DIR"
mkdir -p "$MORPH_OUTPUT_DIR"

echo "Host: $(hostname)"
echo "Project: $PROJECT_DIR"
echo "Candidates: $MORPH_NTRIALS"
echo "Tilts: $MORPH_TILTS"
echo "Boresights: $MORPH_THETA0S"
echo "Output: $MORPH_OUTPUT_DIR"

run_arguments=(
    --ntrials "$MORPH_NTRIALS"
    --output-dir "$MORPH_OUTPUT_DIR"
    --tilts "$MORPH_TILTS"
    --theta0s "$MORPH_THETA0S"
    --julia "$JULIA_BIN"
    --write-batch-events "$MORPH_WRITE_BATCH"
)
case "${MORPH_SKIP_PLOTS,,}" in
    1|true|yes) run_arguments+=(--skip-plots) ;;
    0|false|no) ;;
    *)
        echo "MORPH_SKIP_PLOTS must be true/false (or 1/0); got '$MORPH_SKIP_PLOTS'." >&2
        exit 2
        ;;
esac

# Do not use the local convenience wrapper's default Julia path on OSC.
bash scripts/run_fixed_facet_morphology_scan.sh "${run_arguments[@]}"
