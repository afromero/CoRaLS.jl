#!/usr/bin/env bash
# Source this file on OSC to define a small CoRaLS working environment.
#
# Recommended shared-location usage, with this file placed beside all clones:
#   source /fs/project/PAS0654/machtay_corals/corals_osc_env.sh
#   corals_osc_worktrees
#   corals_osc_activate main
#   corals_osc_activate CoRaLS-osc-slope
#
# A worktree selector can be `main`, a directory name in the shared clone
# directory, or an explicit project path.  The helper also works when sourced
# from `scripts/` inside a checkout.
#
# Explicit-path usage from any directory:
#   source /path/to/CoRaLS.jl/scripts/corals_osc_env.sh
#   corals_osc_activate /path/to/CoRaLS_osc_slope
#
# The functions deliberately do not activate automatically when sourced, so
# merely opening a shell does not change its modules, working directory, or
# Julia depot.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "Source this helper; do not execute it:" >&2
    echo "  source $0 && corals_osc_activate" >&2
    exit 2
fi

_corals_osc_helper_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${_corals_osc_helper_dir}/../Project.toml" ]]; then
    # The helper was sourced from scripts/ inside one checkout.
    _corals_osc_worktree_root="$(cd "${_corals_osc_helper_dir}/../.." && pwd)"
else
    # The helper was placed beside the main clone and its worktrees.
    _corals_osc_worktree_root="${_corals_osc_helper_dir}"
fi

_corals_osc_resolve_worktree() {
    local selector="${1:-}"
    local candidate
    case "$selector" in
        "")
            if [[ -n "${CORALS_PROJECT_DIR:-}" && -f "${CORALS_PROJECT_DIR}/Project.toml" ]]; then
                candidate="$CORALS_PROJECT_DIR"
            else
                candidate="${_corals_osc_worktree_root}/CoRaLS.jl"
            fi
            ;;
        main)
            candidate="${_corals_osc_worktree_root}/CoRaLS.jl"
            ;;
        slope)
            candidate="${_corals_osc_worktree_root}/CoRaLS-osc-slope"
            ;;
        *)
            # Prefer an explicit path; otherwise interpret the selector as a
            # sibling directory name such as CoRaLS-osc-slope.
            if [[ -d "$selector" ]]; then
                candidate="$selector"
            else
                candidate="${_corals_osc_worktree_root}/${selector}"
            fi
            ;;
    esac
    [[ -f "$candidate/Project.toml" ]] || {
        echo "No CoRaLS worktree found for '$selector': $candidate" >&2
        echo "Available worktrees:" >&2
        corals_osc_worktrees >&2
        return 2
    }
    cd "$candidate" && pwd
}

corals_osc_worktrees() {
    local candidate branch
    for candidate in "${_corals_osc_worktree_root}"/*; do
        [[ -f "$candidate/Project.toml" ]] || continue
        branch="$(git -C "$candidate" branch --show-current 2>/dev/null || printf 'not a git checkout')"
        printf '%-32s %s\n' "$(basename "$candidate")" "$branch"
    done
}

corals_osc_activate() {
    local project_dir
    project_dir="$(_corals_osc_resolve_worktree "${1:-}")" || return
    if ! type module >/dev/null 2>&1; then
        echo "The OSC module command is unavailable in this shell." >&2
        return 2
    fi

    module load julia/1.11 python/3.12 || return

    export CORALS_PROJECT_DIR="$project_dir"
    export CORALS_JULIA_DEPOT="${CORALS_JULIA_DEPOT:-${_corals_osc_worktree_root}/julia_depot}"
    mkdir -p "$CORALS_JULIA_DEPOT" || return
    # The final ':' preserves Julia's module-provided system depot while
    # putting registries, packages, artifacts, and compiled caches on /fs.
    export JULIA_DEPOT_PATH="${CORALS_JULIA_DEPOT}:"
    export JULIA_PROJECT="$CORALS_PROJECT_DIR"
    export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-1}"
    # PyCall reads this variable when it is explicitly rebuilt by the setup
    # script.  It must name the module-provided Python, not /usr/bin/python3.
    export PYTHON="${PYTHON:-$(command -v python3)}"

    echo "CoRaLS project: $CORALS_PROJECT_DIR"
    echo "Julia depot:    $CORALS_JULIA_DEPOT"
    echo "Julia:          $(julia --version)"
    echo "Python:         $(python3 --version)"
    echo "PyCall Python:  $PYTHON"
}

corals_osc_cd() {
    [[ -n "${CORALS_PROJECT_DIR:-}" ]] || {
        echo "Run corals_osc_activate first." >&2
        return 2
    }
    cd "$CORALS_PROJECT_DIR"
}

corals_osc_julia() {
    [[ -n "${CORALS_PROJECT_DIR:-}" ]] || {
        echo "Run corals_osc_activate first." >&2
        return 2
    }
    (
        cd "$CORALS_PROJECT_DIR" || exit
        julia --project=. "$@"
    )
}

corals_osc_morphology() {
    [[ -n "${CORALS_PROJECT_DIR:-}" ]] || {
        echo "Run corals_osc_activate first." >&2
        return 2
    }
    (
        cd "$CORALS_PROJECT_DIR" || exit
        bash scripts/run_fixed_facet_morphology_scan.sh --julia "$(command -v julia)" "$@"
    )
}

corals_osc_check_python() {
    python3 -c 'import matplotlib, numpy, pandas; print("Python plotting packages are ready")'
}
