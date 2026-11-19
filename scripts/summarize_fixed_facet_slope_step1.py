#!/usr/bin/env python3
"""Stream a fixed-facet diagnostic CSV into a small binned histogram cache.

The raw paired-event table intentionally repeats payload/geometry columns for
every route, boresight, and facet tilt.  This script reads it in chunks and
writes only the histogram counts needed by the distribution plotter.  It never
loads the full raw table into memory.
"""

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


FIELDS = [
    ("raw_exit_angle_theta_r_deg", "theta_r", (0, 95)),
    ("theta_i_deg", "theta_i", (0, 50)),
    ("theta_emit_deg", "theta_emit", (0, 50)),
    ("psi_deg", "psi", (0, 180)),
    ("payload_elevation_signal_deg", "payload_elevation", (-35, 0)),
    ("payload_zenith_at_surface_deg", "payload_zenith", (0, 95)),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/terrain_diagnostics/fixed_facet_slope_step1.csv")
    parser.add_argument("--output", default="data/terrain_diagnostics/fixed_facet_slope_step1_histograms.csv")
    parser.add_argument("--pairing-output",
                        help="small exact gained/lost summary; default derives from --output")
    parser.add_argument("--chunksize", type=int, default=200_000)
    args = parser.parse_args()
    if args.chunksize <= 0:
        raise ValueError("--chunksize must be positive")

    usecols = ["event_id", "route", "facet_tilt_deg", "lpda_theta0_deg", "status", "triggered"] + [f[0] for f in FIELDS]
    counts = defaultdict(lambda: np.zeros(60, dtype=np.int64))
    pairing = defaultdict(lambda: np.zeros(5, dtype=np.int64))

    def accumulate_pairing(frame):
        """Accumulate matched-event counts from complete event-ID groups only."""
        if frame.empty:
            return
        for (route, theta0), panel in frame.groupby(["route", "lpda_theta0_deg"], sort=False):
            for selection, passed in (
                ("valid", panel.status.eq("signal")),
                ("triggered", panel.status.eq("signal") & panel.triggered.astype(str).str.lower().eq("true")),
            ):
                table = panel.loc[:, ["event_id", "facet_tilt_deg"]].copy()
                table["passed"] = passed.to_numpy()
                baseline = table[np.isclose(table.facet_tilt_deg, 0.0)].set_index("event_id").passed
                if baseline.empty:
                    continue
                for tilt, current_table in table.groupby("facet_tilt_deg", sort=False):
                    current = current_table.set_index("event_id").passed
                    shared = baseline.index.intersection(current.index)
                    if shared.empty:
                        continue
                    base = baseline.loc[shared].to_numpy(dtype=bool)
                    value = current.loc[shared].to_numpy(dtype=bool)
                    key = (str(route), float(theta0), selection, float(tilt))
                    pairing[key] += np.array([
                        base.sum(),
                        np.logical_and(base, value).sum(),
                        np.logical_and(base, ~value).sum(),
                        np.logical_and(~base, value).sum(),
                        value.sum(),
                    ])

    nrows = 0
    carry = pd.DataFrame()
    previous_last_event_id = None
    for chunk in pd.read_csv(args.input, usecols=usecols, chunksize=args.chunksize):
        # Exact pairing below relies on the recorder's event-major output:
        # every route/tilt/boresight row for one event is consecutive.  Fail
        # explicitly rather than silently making an incomplete pairing plot if
        # a hand-edited or externally sorted table violates that contract.
        event_ids = chunk.event_id.to_numpy()
        if np.any(event_ids[1:] < event_ids[:-1]) or (
            previous_last_event_id is not None and event_ids[0] < previous_last_event_id
        ):
            raise ValueError(
                "event_id must be nondecreasing: pairing requires event-major rows "
                "from record_fixed_facet_slope_step1.jl"
            )
        previous_last_event_id = event_ids[-1]
        nrows += len(chunk)
        valid = chunk.status.eq("signal")
        triggered = valid & chunk.triggered.astype(str).str.lower().eq("true")
        for field, token, value_range in FIELDS:
            values = pd.to_numeric(chunk[field], errors="coerce")
            for selection, mask in (("valid", valid), ("triggered", triggered)):
                selected = chunk.loc[mask, ["route", "facet_tilt_deg", "lpda_theta0_deg"]].copy()
                selected["value"] = values[mask]
                selected = selected.dropna(subset=["value"])
                selected = selected[selected.value.between(*value_range)]
                for (route, tilt, theta0), group in selected.groupby(
                    ["route", "facet_tilt_deg", "lpda_theta0_deg"], sort=False
                ):
                    key = (str(route), float(tilt), float(theta0), selection, token)
                    counts[key] += np.histogram(group.value.to_numpy(), bins=60, range=value_range)[0]
        # The recorder writes all rows for an event ID contiguously.  Preserve
        # the final (possibly split) event as carry; every earlier ID can be
        # compared against its own zero-tilt row exactly.
        pairing_chunk = chunk[["event_id", "route", "facet_tilt_deg", "lpda_theta0_deg", "status", "triggered"]]
        if not carry.empty:
            pairing_chunk = pd.concat((carry, pairing_chunk), ignore_index=True)
        last_event_id = pairing_chunk.event_id.iloc[-1]
        complete = pairing_chunk[pairing_chunk.event_id != last_event_id]
        carry = pairing_chunk[pairing_chunk.event_id == last_event_id].copy()
        accumulate_pairing(complete)
        print(f"Read {nrows:,} rows", end="\r", flush=True)

    accumulate_pairing(carry)

    rows = []
    for (route, tilt, theta0, selection, token), count in sorted(counts.items()):
        value_range = next(item[2] for item in FIELDS if item[1] == token)
        edges = np.linspace(*value_range, len(count) + 1)
        rows.extend(
            {
                "route": route,
                "facet_tilt_deg": tilt,
                "lpda_theta0_deg": theta0,
                "selection": selection,
                "field": token,
                "bin_left_deg": left,
                "bin_right_deg": right,
                "count": int(value),
            }
            for left, right, value in zip(edges[:-1], edges[1:], count)
        )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False)
    pairing_output = Path(args.pairing_output) if args.pairing_output else output.with_name(
        output.stem + "_pairing_counts.csv"
    )
    pairing_rows = [
        {
            "route": route,
            "lpda_theta0_deg": theta0,
            "selection": selection,
            "facet_tilt_deg": tilt,
            "baseline_count": int(value[0]),
            "retained_count": int(value[1]),
            "lost_count": int(value[2]),
            "gained_count": int(value[3]),
            "total_count": int(value[4]),
        }
        for (route, theta0, selection, tilt), value in sorted(pairing.items())
    ]
    pairing_output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(pairing_rows).to_csv(pairing_output, index=False)
    print(f"\nWrote {output} ({len(rows):,} histogram bins from {nrows:,} raw rows)")
    print(f"Wrote {pairing_output} ({len(pairing_rows):,} paired-count rows)")


if __name__ == "__main__":
    main()
