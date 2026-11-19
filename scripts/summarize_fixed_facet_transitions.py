#!/usr/bin/env python3
"""Stream matched fixed-facet events into adjacent-tilt gained/lost summaries.

For a lost event, physical-angle histograms use the last valid (lower-tilt)
row.  For a gained event they use the newly valid (higher-tilt) row.  The
failure-status table separately records why the opposite endpoint failed.

The raw CSV is event-major: all route/tilt rows for one event ID are
contiguous.  This script preserves only binned distributions, transition
counts/failure reasons, and a small deterministic example set.
"""

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


FIELDS = [
    ("raw_exit_angle_theta_r_deg", "theta_r", (0.0, 95.0)),
    ("theta_i_deg", "theta_i", (0.0, 50.0)),
    ("theta_emit_deg", "theta_emit", (0.0, 50.0)),
    ("theta_ice_deg", "theta_ice", (0.0, 60.0)),
    ("psi_deg", "psi", (0.0, 180.0)),
    ("payload_elevation_signal_deg", "payload_elevation", (-90.0, 0.0)),
    ("payload_zenith_at_surface_deg", "payload_zenith", (0.0, 95.0)),
]
EXAMPLE_COLUMNS = [
    "event_id", "route", "lpda_theta0_deg", "selection", "change",
    "tilt_from_deg", "tilt_to_deg", "facet_aspect_deg",
    "surface_lat_deg", "surface_lon_deg", "payload_lat_deg",
    "payload_lon_deg", "payload_altitude_km", "previous_status",
    "next_status", "previous_triggered", "next_triggered",
    "previous_raw_exit_angle_theta_r_deg", "next_raw_exit_angle_theta_r_deg",
    "previous_theta_i_deg", "next_theta_i_deg",
    "previous_theta_emit_deg", "next_theta_emit_deg",
    "previous_theta_ice_deg", "next_theta_ice_deg",
    "previous_psi_deg", "next_psi_deg",
    "previous_payload_elevation_signal_deg", "next_payload_elevation_signal_deg",
]


def bool_triggered(series):
    return series.astype(str).str.lower().eq("true")


def score(event_id, key):
    """Stable 64-bit mixing for deterministic representative-event sampling."""
    key_value = 0
    for byte in repr(key).encode("utf-8"):
        key_value = (key_value * 131 + byte) & ((1 << 64) - 1)
    value = (int(event_id) + 0x9E3779B97F4A7C15 + key_value) & ((1 << 64) - 1)
    value ^= value >> 30
    value = (value * 0xBF58476D1CE4E5B9) & ((1 << 64) - 1)
    value ^= value >> 27
    value = (value * 0x94D049BB133111EB) & ((1 << 64) - 1)
    return value ^ (value >> 31)


def endpoint_record(event_id, route, theta0, selection, change, previous, current):
    """Keep both endpoint states; physical histograms choose one by change."""
    record = {
        "event_id": int(event_id), "route": str(route), "lpda_theta0_deg": float(theta0),
        "selection": selection, "change": change,
        "tilt_from_deg": float(previous.facet_tilt_deg),
        "tilt_to_deg": float(current.facet_tilt_deg),
        "facet_aspect_deg": previous.facet_aspect_deg,
        "surface_lat_deg": previous.surface_lat_deg,
        "surface_lon_deg": previous.surface_lon_deg,
        "payload_lat_deg": previous.payload_lat_deg,
        "payload_lon_deg": previous.payload_lon_deg,
        "payload_altitude_km": previous.payload_altitude_km,
        "previous_status": previous.status, "next_status": current.status,
        "previous_triggered": bool(previous.triggered_bool),
        "next_triggered": bool(current.triggered_bool),
    }
    for field, _, _ in FIELDS:
        record[f"previous_{field}"] = getattr(previous, field)
        record[f"next_{field}"] = getattr(current, field)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/terrain_diagnostics/fixed_facet_slope_step1.csv")
    parser.add_argument("--output-prefix",
                        default="data/terrain_diagnostics/fixed_facet_slope_step1_transitions")
    parser.add_argument("--chunksize", type=int, default=200_000)
    parser.add_argument("--examples-per-category", type=int, default=24)
    args = parser.parse_args()
    if args.chunksize <= 0 or args.examples_per_category <= 0:
        raise ValueError("--chunksize and --examples-per-category must be positive")

    usecols = [
        "event_id", "route", "facet_tilt_deg", "facet_aspect_deg", "lpda_theta0_deg",
        "status", "triggered", "surface_lat_deg", "surface_lon_deg", "payload_lat_deg",
        "payload_lon_deg", "payload_altitude_km",
    ] + [field for field, _, _ in FIELDS]
    counts = defaultdict(lambda: np.zeros(60, dtype=np.int64))
    transitions = defaultdict(int)
    failures = defaultdict(int)
    examples = defaultdict(list)
    ranges = {token: value_range for _, token, value_range in FIELDS}
    carry = pd.DataFrame()
    previous_last_event_id = None
    nrows = 0

    def update_examples(key, records):
        """Keep the first event-major rows; CR/event draws are already randomized."""
        remaining = args.examples_per_category - len(examples[key])
        if remaining > 0:
            examples[key].extend(records[:remaining])

    def process_events(frame):
        if frame.empty:
            return
        frame = frame.copy()
        frame["triggered_bool"] = bool_triggered(frame.triggered)
        tilts = sorted(frame.facet_tilt_deg.unique())
        join_columns = ["event_id", "route", "lpda_theta0_deg"]
        for lo, hi in zip(tilts, tilts[1:]):
            previous = frame[np.isclose(frame.facet_tilt_deg, lo)]
            current = frame[np.isclose(frame.facet_tilt_deg, hi)]
            paired = previous.merge(current, on=join_columns, suffixes=("_previous", "_next"),
                                    validate="one_to_one")
            for selection, prev_pass, curr_pass in (
                ("valid", paired.status_previous.eq("signal"), paired.status_next.eq("signal")),
                ("triggered", paired.status_previous.eq("signal") & paired.triggered_bool_previous,
                 paired.status_next.eq("signal") & paired.triggered_bool_next),
            ):
                for change, mask, endpoint, opposite in (
                    ("gained", curr_pass & ~prev_pass, "next", "previous"),
                    ("lost", prev_pass & ~curr_pass, "previous", "next"),
                ):
                    changed = paired.loc[mask]
                    for (route, theta0), group in changed.groupby(["route", "lpda_theta0_deg"], sort=False):
                        key = (str(route), float(theta0), selection, float(lo), float(hi), change)
                        transitions[key] += len(group)
                        for status, amount in group[f"status_{opposite}"].value_counts().items():
                            failures[(key, str(status))] += int(amount)
                        for field, token, value_range in FIELDS:
                            values = pd.to_numeric(group[f"{field}_{endpoint}"], errors="coerce").dropna()
                            values = values[values.between(*value_range)]
                            if not values.empty:
                                counts[(key, token)] += np.histogram(values.to_numpy(), bins=60,
                                                                      range=value_range)[0]
                        records = []
                        for row in group.head(args.examples_per_category).itertuples(index=False):
                            record = {
                                "event_id": int(row.event_id), "route": str(route),
                                "lpda_theta0_deg": float(theta0), "selection": selection,
                                "change": change, "tilt_from_deg": float(lo), "tilt_to_deg": float(hi),
                                "facet_aspect_deg": row.facet_aspect_deg_previous,
                                "surface_lat_deg": row.surface_lat_deg_previous,
                                "surface_lon_deg": row.surface_lon_deg_previous,
                                "payload_lat_deg": row.payload_lat_deg_previous,
                                "payload_lon_deg": row.payload_lon_deg_previous,
                                "payload_altitude_km": row.payload_altitude_km_previous,
                                "previous_status": row.status_previous, "next_status": row.status_next,
                                "previous_triggered": bool(row.triggered_bool_previous),
                                "next_triggered": bool(row.triggered_bool_next),
                            }
                            for field, _, _ in FIELDS:
                                record[f"previous_{field}"] = getattr(row, f"{field}_previous")
                                record[f"next_{field}"] = getattr(row, f"{field}_next")
                            records.append(record)
                        update_examples(key, records)

    for chunk in pd.read_csv(args.input, usecols=usecols, chunksize=args.chunksize):
        event_ids = chunk.event_id.to_numpy()
        if np.any(event_ids[1:] < event_ids[:-1]) or (
            previous_last_event_id is not None and event_ids[0] < previous_last_event_id
        ):
            raise ValueError("event_id must be nondecreasing and event-major")
        previous_last_event_id = event_ids[-1]
        nrows += len(chunk)
        if not carry.empty:
            chunk = pd.concat((carry, chunk), ignore_index=True)
        last_event_id = chunk.event_id.iloc[-1]
        complete = chunk[chunk.event_id != last_event_id]
        carry = chunk[chunk.event_id == last_event_id].copy()
        process_events(complete)
        print(f"Read {nrows:,} rows", end="\r", flush=True)
    process_events(carry)

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    histogram_rows = []
    for (key, token), count in sorted(counts.items()):
        route, theta0, selection, lo, hi, change = key
        edges = np.linspace(*ranges[token], len(count) + 1)
        for left, right, value in zip(edges[:-1], edges[1:], count):
            histogram_rows.append({
                "route": route, "lpda_theta0_deg": theta0, "selection": selection,
                "tilt_from_deg": lo, "tilt_to_deg": hi, "change": change, "field": token,
                "bin_left_deg": left, "bin_right_deg": right, "count": int(value),
            })
    hist_path = prefix.with_name(prefix.name + "_histograms.csv")
    pd.DataFrame(histogram_rows).to_csv(hist_path, index=False)

    transition_rows = []
    for key, count in sorted(transitions.items()):
        route, theta0, selection, lo, hi, change = key
        transition_rows.append({
            "route": route, "lpda_theta0_deg": theta0, "selection": selection,
            "tilt_from_deg": lo, "tilt_to_deg": hi, "change": change, "count": count,
        })
    transition_path = prefix.with_name(prefix.name + "_counts.csv")
    pd.DataFrame(transition_rows).to_csv(transition_path, index=False)

    failure_rows = []
    for (key, status), count in sorted(failures.items()):
        route, theta0, selection, lo, hi, change = key
        failure_rows.append({
            "route": route, "lpda_theta0_deg": theta0, "selection": selection,
            "tilt_from_deg": lo, "tilt_to_deg": hi, "change": change,
            "opposite_endpoint_status": status, "count": count,
        })
    failure_path = prefix.with_name(prefix.name + "_failure_statuses.csv")
    pd.DataFrame(failure_rows).to_csv(failure_path, index=False)

    example_rows = []
    for records in examples.values():
        example_rows.extend(records)
    example_path = prefix.with_name(prefix.name + "_examples.csv")
    pd.DataFrame(example_rows, columns=EXAMPLE_COLUMNS).sort_values(
        ["route", "lpda_theta0_deg", "selection", "tilt_from_deg", "change", "event_id"]
    ).to_csv(example_path, index=False)
    print(f"\nWrote {hist_path} ({len(histogram_rows):,} bins)")
    print(f"Wrote {transition_path} ({len(transition_rows):,} transition rows)")
    print(f"Wrote {failure_path} ({len(failure_rows):,} failure-status rows)")
    print(f"Wrote {example_path} ({len(example_rows):,} representative event rows)")


if __name__ == "__main__":
    main()
