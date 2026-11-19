#!/usr/bin/env bash
# Run the paired fixed-facet morphology diagnostic and its compact analyses.
#
# The simulation shares each candidate event geometry across every facet tilt,
# so differences between tilts can be interpreted event by event.  The raw CSV
# is deliberately detailed and can be several GB; the subsequent Python steps
# turn it into small histogram/count caches and PNG figures without loading the
# full CSV into memory.
#
# Default scan: 500,000 shared candidates, tilts 0--10 degrees in one-degree
# steps, and the nominal -48 degree LPDA boresight.  This produces about
# 4.7 GB of raw CSV, plus only a few MB of summarized products.
#
# Usage:
#   bash scripts/run_fixed_facet_morphology_scan.sh [options]
#
# Options:
#   --ntrials N        Shared candidate geometries (default: 500000)
#   --output-dir PATH  Isolated directory for all products
#   --tilts LIST       Comma-separated facet tilts in degrees
#   --theta0s LIST     Comma-separated LPDA boresights in degrees
#   --julia PATH       Julia executable
#   --write-batch-events N
#                       Candidate events retained before appending CSV rows
#   --skip-plots        Build raw data and compact summaries, but do not invoke
#                       matplotlib plotting scripts
#   --skip-simulation  Rebuild summaries/figures from an existing events.csv
#   -h, --help         Show this message
#
# Examples:
#   # Default 0--10 degree, one-degree scan.
#   bash scripts/run_fixed_facet_morphology_scan.sh
#
#   # Quicker exploratory scan, written somewhere else.
#   bash scripts/run_fixed_facet_morphology_scan.sh \
#     --ntrials 100000 --output-dir data/terrain_diagnostics/fixed_facet_quick
#
#   # Regenerate compact products without rerunning the Monte Carlo.
#   bash scripts/run_fixed_facet_morphology_scan.sh --skip-simulation

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

ntrials=500000
output_dir="data/terrain_diagnostics/fixed_facet_0to10_1deg"
tilts="0,1,2,3,4,5,6,7,8,9,10"
theta0s="-48"
julia_bin="${JULIA_BIN:-/home/machtay/julia-1.11.9/bin/julia}"
write_batch_events=5000
skip_plots=false
skip_simulation=false

usage() {
    sed -n '2,/^set -euo pipefail$/p' "$0" | sed '$d; s/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ntrials)
            ntrials="$2"; shift 2 ;;
        --output-dir)
            output_dir="$2"; shift 2 ;;
        --tilts)
            tilts="$2"; shift 2 ;;
        --theta0s)
            theta0s="$2"; shift 2 ;;
        --julia)
            julia_bin="$2"; shift 2 ;;
        --write-batch-events)
            write_batch_events="$2"; shift 2 ;;
        --skip-plots)
            skip_plots=true; shift ;;
        --skip-simulation)
            skip_simulation=true; shift ;;
        -h|--help)
            usage; exit 0 ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2 ;;
    esac
done

if [[ ! "$ntrials" =~ ^[1-9][0-9]*$ ]]; then
    echo "--ntrials must be a positive integer; got '$ntrials'" >&2
    exit 2
fi
if [[ ! "$write_batch_events" =~ ^[1-9][0-9]*$ ]]; then
    echo "--write-batch-events must be a positive integer; got '$write_batch_events'" >&2
    exit 2
fi
if [[ ! -x "$julia_bin" ]]; then
    echo "Julia executable not found: $julia_bin" >&2
    echo "Pass --julia PATH or set JULIA_BIN=PATH." >&2
    exit 2
fi

events="$output_dir/events.csv"
histograms="$output_dir/histograms.csv"
pairing="$output_dir/histograms_pairing_counts.csv"
transitions_prefix="$output_dir/transitions"
figures="$output_dir/figures"
mkdir -p "$figures"

echo "Fixed-facet morphology scan"
echo "  candidates: $ntrials"
echo "  tilts:      $tilts deg"
echo "  boresights: $theta0s deg"
echo "  output:     $output_dir"

if [[ "$skip_simulation" == false ]]; then
    if [[ -e "$events" ]]; then
        echo "Refusing to overwrite existing raw event table: $events" >&2
        echo "Choose another --output-dir, move the file, or use --skip-simulation." >&2
        exit 1
    fi
    "$julia_bin" --project=. scripts/record_fixed_facet_slope_step1.jl \
        --ntrials "$ntrials" \
        --facet-tilts "$tilts" \
        --theta0s "$theta0s" \
        --write-batch-events "$write_batch_events" \
        --output "$events"
elif [[ ! -f "$events" ]]; then
    echo "--skip-simulation requires an existing event table: $events" >&2
    exit 1
fi

python3 scripts/summarize_fixed_facet_slope_step1.py \
    --input "$events" \
    --output "$histograms"

python3 scripts/summarize_fixed_facet_transitions.py \
    --input "$events" \
    --output-prefix "$transitions_prefix"

if [[ "$skip_plots" == false ]]; then
    python3 scripts/plot_fixed_facet_slope_step1.py \
        --histogram-input "$histograms" \
        --pairing-summary "$pairing" \
        --output-dir "$figures" \
        --pairing-output "$figures/pairing.png" \
        --layout grid \
        --grid-cols 3

    # Render both paths.  Each route has valid-scalar and LPDA-triggered grids for
    # every adjacent tilt pair, with lost events in red and gained events in green.
    for route in direct reflected; do
        python3 scripts/plot_fixed_facet_transitions.py \
            --input "${transitions_prefix}_histograms.csv" \
            --output-dir "$figures/transitions" \
            --route "$route" \
            --normalization density
    done
else
    echo "Skipping matplotlib figure generation; raw data and compact summaries were retained."
fi

echo
echo "Complete.  Raw events: $events"
echo "Compact all-tilt cache: $histograms"
echo "Compact adjacent-transition cache: ${transitions_prefix}_histograms.csv"
echo "Figures: $figures"
