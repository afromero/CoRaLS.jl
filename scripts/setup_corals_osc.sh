#!/usr/bin/env bash
# One-time (or safe-to-repeat) CoRaLS setup on OSC.
#
# This script loads the standard Julia/Python modules, establishes a shared
# Julia depot beside the selected worktree, instantiates and precompiles the
# project, checks the Python plotting dependencies, and creates Slurm logs.
#
# Usage:
#   bash scripts/setup_corals_osc.sh
#   bash scripts/setup_corals_osc.sh --project-dir /fs/project/.../CoRaLS-osc-slope
#   bash scripts/setup_corals_osc.sh --depot /fs/project/.../julia_depot

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "$script_dir/.." && pwd)"
depot_dir=""

usage() {
    sed -n '2,/^set -euo pipefail$/p' "$0" | sed '$d; s/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --project-dir)
            project_dir="$2"; shift 2 ;;
        --depot)
            depot_dir="$2"; shift 2 ;;
        -h|--help)
            usage; exit 0 ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2 ;;
    esac
done

project_dir="$(cd "$project_dir" && pwd)"
[[ -f "$project_dir/Project.toml" ]] || {
    echo "Not a CoRaLS project: $project_dir" >&2
    exit 2
}

if [[ -n "$depot_dir" ]]; then
    export CORALS_JULIA_DEPOT="$depot_dir"
fi
export CORALS_PROJECT_DIR="$project_dir"

# shellcheck source=corals_osc_env.sh
source "$script_dir/corals_osc_env.sh"
corals_osc_activate "$project_dir"
corals_osc_check_python

mkdir -p "$CORALS_PROJECT_DIR/slurm/out" "$CORALS_PROJECT_DIR/slurm/err"

echo "Instantiating and precompiling Julia dependencies..."
corals_osc_julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); Pkg.status()'

cat <<EOF

OSC setup is complete.

For later interactive sessions:
  source $CORALS_PROJECT_DIR/scripts/corals_osc_env.sh
  corals_osc_activate
  corals_osc_cd

Useful commands after activation:
  corals_osc_julia -e 'using CoRaLS; println("CoRaLS ready")'
  corals_osc_morphology --ntrials 1000 --output-dir /fs/project/.../quick_check
EOF
