#!/usr/bin/env python3
"""Plot paired fixed-facet direct/R distributions and event membership changes."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_DIR = Path("data/terrain_diagnostics")
FIELDS = [
    ("raw_exit_angle_theta_r_deg", "theta_r", (0, 95), r"Local exit angle $\theta_r$ [deg]"),
    ("theta_i_deg", "theta_i", (0, 50), r"Surface-side $\theta_i$ [deg]"),
    ("theta_emit_deg", "theta_emit", (0, 50), r"Ray angle at scalar source [deg]"),
    ("psi_deg", "psi", (0, 180), r"Shower off-axis angle $\psi$ [deg]"),
    ("payload_elevation_signal_deg", "payload_elevation", (-35, 0),
     r"Payload elevation $\theta_{\rm el}$ [deg]"),
    ("payload_zenith_at_surface_deg", "payload_zenith", (0, 95), "Payload zenith at surface [deg]"),
]


def count_hist(axis, values, value_range, color, label, style, *, normalization="counts", fill=False):
    values = values.dropna()
    values = values[(values >= value_range[0]) & (values <= value_range[1])]
    if values.empty:
        return
    edges = np.linspace(*value_range, 61)
    count, _ = np.histogram(values, bins=edges)
    if normalization == "density":
        integral = count.sum() * np.diff(edges)[0]
        if integral <= 0:
            return
        count = count / integral
    elif normalization != "counts":
        raise ValueError(f"unknown normalization: {normalization}")
    axis.stairs(count, edges, color=color, linestyle=style, linewidth=1.35, alpha=0.86, label=label)
    if fill:
        axis.stairs(count, edges, color=color, fill=True, alpha=0.10)


def selection_title(selection):
    return "valid scalar solutions (boresight-independent)" if selection == "valid" else "LPDA-triggered solutions"


def plot_distribution_axis(axis, panel, selection, field, value_range, xlabel, field_label, theta0=None,
                           normalization="counts"):
    """Populate one axis with the matched-facet histogram for a single quantity."""
    tilts = sorted(panel.facet_tilt_deg.unique())
    colors = plt.cm.viridis_r(np.linspace(0.08, 0.92, len(tilts)))
    for tilt, color in zip(tilts, colors):
        group = panel[np.isclose(panel.facet_tilt_deg, tilt)]
        selected = group[group.status == "signal"]
        if selection == "triggered":
            selected = selected[selected.triggered]
        count_hist(axis, selected[field], value_range, color, f"{tilt:g}°", "-",
                   normalization=normalization, fill=selection == "triggered")
    title = field_label
    if theta0 is not None:
        title += fr" · $\theta_0={theta0:g}^\circ$"
    ylabel = "Events per angle bin" if normalization == "counts" else "Normalized density [deg$^{-1}$]"
    axis.set(xlabel=xlabel, ylabel=ylabel, title=title)
    axis.grid(alpha=0.25)
    if field == "raw_exit_angle_theta_r_deg":
        axis.axvline(90, color="0.2", linestyle=":", linewidth=1)


def add_legend(axis):
    handles, labels = axis.get_legend_handles_labels()
    if handles:
        axis.legend(handles, labels, title="Facet tilt", fontsize=8, title_fontsize=9, ncol=2, loc="best")


def distribution_figure(panel, output, route, selection, field, field_label, value_range, xlabel, theta0=None,
                        normalization="counts"):
    """Plot one angle quantity, route, selection, and (when relevant) boresight."""
    figure, axis = plt.subplots(figsize=(7.8, 5.2))
    plot_distribution_axis(axis, panel, selection, field, value_range, xlabel, field_label, theta0, normalization)
    add_legend(axis)
    norm_title = "Raw event counts" if normalization == "counts" else "Each curve normalized to unit integral"
    figure.suptitle(norm_title + ": paired fixed-facet slopes reuse the same event IDs, payload positions,\n"
                     "CR axes, and facet azimuths across tilt magnitudes.", fontsize=11)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)
    print(f"Wrote {output}")


def distribution_grid_figure(panel, output, route, selection, ncols, theta0=None, normalization="counts"):
    """Plot all angle quantities in a compact multi-row figure."""
    nrows = int(np.ceil(len(FIELDS) / ncols))
    figure, axes = plt.subplots(nrows, ncols, figsize=(7.2 * ncols, 4.6 * nrows), squeeze=False)
    flat_axes = axes.ravel()
    for axis, (field, field_token, value_range, xlabel) in zip(flat_axes, FIELDS):
        plot_distribution_axis(axis, panel, selection, field, value_range, xlabel, field_token, theta0,
                               normalization)
    for axis in flat_axes[len(FIELDS):]:
        axis.set_visible(False)
    add_legend(axes[0, 0])
    title = f"{route} · {selection_title(selection)}"
    if theta0 is not None:
        title += fr" · $\theta_0={theta0:g}^\circ$"
    norm_title = "Raw event counts" if normalization == "counts" else "Each curve normalized to unit integral"
    figure.suptitle(title + " · " + norm_title +
                     "\nMatched event IDs, payload positions, CR axes, and facet azimuths across all tilts.",
                     fontsize=13)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)
    print(f"Wrote {output}")


def plot_histogram_grid(cache, output, route, selection, ncols, theta0=None, normalization="counts"):
    """Render a grid from the compact histogram cache made by the summarizer."""
    nrows = int(np.ceil(len(FIELDS) / ncols))
    figure, axes = plt.subplots(nrows, ncols, figsize=(7.2 * ncols, 4.6 * nrows), squeeze=False)
    flat_axes = axes.ravel()
    tilts = sorted(cache.facet_tilt_deg.unique())
    colors = plt.cm.viridis_r(np.linspace(0.08, 0.92, len(tilts)))
    for axis, (_, token, value_range, xlabel) in zip(flat_axes, FIELDS):
        panel = cache[(cache.route == route) & (cache.selection == selection) & (cache.field == token)]
        if theta0 is not None:
            panel = panel[np.isclose(panel.lpda_theta0_deg, theta0)]
        for tilt, color in zip(tilts, colors):
            group = panel[np.isclose(panel.facet_tilt_deg, tilt)].sort_values("bin_left_deg")
            if group.empty:
                continue
            count = group["count"].to_numpy(dtype=float)
            edges = np.r_[group.bin_left_deg.to_numpy(), group.bin_right_deg.iloc[-1]]
            if normalization == "density":
                integral = count.sum() * np.diff(edges)[0]
                if integral > 0:
                    count /= integral
            axis.stairs(count, edges, color=color, linewidth=1.35, alpha=0.86, label=f"{tilt:g}°")
            if selection == "triggered":
                axis.stairs(count, edges, color=color, fill=True, alpha=0.10)
        ylabel = "Events per angle bin" if normalization == "counts" else "Normalized density [deg$^{-1}$]"
        title = token
        if theta0 is not None:
            title += fr" · $\theta_0={theta0:g}^\circ$"
        axis.set(xlabel=xlabel, ylabel=ylabel, title=title, xlim=value_range)
        axis.grid(alpha=0.25)
        if token == "theta_r":
            axis.axvline(90, color="0.2", linestyle=":", linewidth=1)
    for axis in flat_axes[len(FIELDS):]:
        axis.set_visible(False)
    add_legend(axes[0, 0])
    title = f"{route} · {selection_title(selection)}"
    if theta0 is not None:
        title += fr" · $\theta_0={theta0:g}^\circ$"
    norm_title = "Raw event counts" if normalization == "counts" else "Each curve normalized to unit integral"
    figure.suptitle(title + " · " + norm_title + "\nStreaming histogram cache; no raw event table loaded.", fontsize=13)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)
    print(f"Wrote {output}")


def plot_from_histogram_cache(cache, output_dir, layout, ncols, normalizations):
    routes = [route for route in ("direct", "reflected") if route in set(cache.route)]
    for route in routes:
        route_cache = cache[cache.route == route]
        reference_theta0 = sorted(route_cache.lpda_theta0_deg.unique())[0]
        for normalization in normalizations:
            suffix = "" if normalization == "counts" else "_density"
            # The valid scalar distributions are repeated for every boresight;
            # select one identical copy.
            valid = route_cache[np.isclose(route_cache.lpda_theta0_deg, reference_theta0)]
            if layout != "grid":
                raise ValueError("--histogram-input currently supports --layout grid only")
            plot_histogram_grid(valid, output_dir / f"fixed_facet_slope_step1_{route}_valid{suffix}.png",
                                route, "valid", ncols, normalization=normalization)
            for theta0 in sorted(route_cache.lpda_theta0_deg.unique()):
                theta0_token = f"theta0_{theta0:+g}".replace("+", "p").replace("-", "m")
                plot_histogram_grid(route_cache,
                                    output_dir / f"fixed_facet_slope_step1_{route}_triggered_{theta0_token}{suffix}.png",
                                    route, "triggered", ncols, theta0=theta0, normalization=normalization)


def membership(frame, output):
    """Compare valid-path and trigger membership with each event's zero-tilt result."""
    tilts = sorted(frame.facet_tilt_deg.unique())
    if not any(np.isclose(tilts, 0.0)):
        raise ValueError("the paired-change plot requires facet tilt 0")
    boresights = sorted(frame.lpda_theta0_deg.unique())
    routes = [route for route in ("direct", "reflected") if route in set(frame.route)]
    outcomes = [
        ("valid scalar path", "valid", lambda panel: panel.status.eq("signal")),
        ("LPDA triggered", "triggered", lambda panel: panel.triggered.astype(bool)),
    ]
    for theta0 in boresights:
        figure, axes = plt.subplots(len(routes), len(outcomes),
                                    figsize=(10.6, 7.4), squeeze=False, sharex=True)
        for route_index, route in enumerate(routes):
            panel = frame[(frame.route == route) & np.isclose(frame.lpda_theta0_deg, theta0)]
            for outcome_index, (outcome_title, column, selector) in enumerate(outcomes):
                axis = axes[route_index, outcome_index]
                membership_frame = panel.assign(**{column: selector(panel)})
                pivot = membership_frame.pivot(index="event_id", columns="facet_tilt_deg", values=column)
                baseline = pivot[0.0].fillna(False).astype(bool)
                lost, gained, totals = [], [], []
                for tilt in tilts:
                    selected = pivot[tilt].fillna(False).astype(bool)
                    lost.append((baseline & ~selected).sum())
                    gained.append((~baseline & selected).sum())
                    totals.append(selected.sum())
                axis.axhline(baseline.sum(), color="0.35", linestyle="--", linewidth=1.2,
                             label="flat-facet count")
                axis.bar(tilts, gained, width=2.7, color="#2ca25f", alpha=0.8, label="gained vs flat")
                axis.bar(tilts, -np.asarray(lost), width=2.7, color="#de2d26", alpha=0.8, label="lost vs flat")
                axis.plot(tilts, totals, "ko-", linewidth=1.4, markersize=3.5, label="total count")
                axis.set(title=f"{route} · {outcome_title} · $\\theta_0={theta0:g}^\\circ$",
                         xlabel="Fixed facet tilt [deg]", ylabel="Paired event count")
                axis.grid(alpha=0.25)
                axis.text(0.03, 0.04, "Green/red bars classify event IDs\nrelative to the same flat-facet draw.",
                          transform=axis.transAxes, fontsize=7.2, va="bottom")
        axes[0, 0].legend(fontsize=7, loc="best")
        figure.suptitle(f"Matched-event changes as fixed facet tilt increases · $\\theta_0={theta0:g}^\\circ$",
                         fontsize=13)
        figure.tight_layout()
        figure_output = pairing_output_for_boresight(output, theta0, len(boresights) > 1)
        figure_output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(figure_output, dpi=180)
        plt.close(figure)
        print(f"Wrote {figure_output}")


def pairing_output_for_boresight(output, theta0, multiple_boresights):
    """Avoid a tall multi-boresight panel: use one readable 2x2 figure each."""
    if not multiple_boresights:
        return output
    token = f"theta0_{theta0:+g}".replace("+", "p").replace("-", "m")
    return output.with_name(f"{output.stem}_{token}{output.suffix}")


def membership_from_summary(summary, output):
    """Render the pairing figure from the tiny exact streaming-count summary."""
    boresights = sorted(summary.lpda_theta0_deg.unique())
    routes = [route for route in ("direct", "reflected") if route in set(summary.route)]
    outcomes = [("valid scalar path", "valid"), ("LPDA triggered", "triggered")]
    for theta0 in boresights:
        figure, axes = plt.subplots(len(routes), len(outcomes),
                                    figsize=(10.6, 7.4), squeeze=False, sharex=True)
        for route_index, route in enumerate(routes):
            panel = summary[(summary.route == route) & np.isclose(summary.lpda_theta0_deg, theta0)]
            for outcome_index, (outcome_title, selection) in enumerate(outcomes):
                axis = axes[route_index, outcome_index]
                rows = panel[panel.selection == selection].sort_values("facet_tilt_deg")
                if rows.empty:
                    continue
                tilts = rows.facet_tilt_deg.to_numpy()
                baseline = rows.baseline_count.iloc[0]
                axis.axhline(baseline, color="0.35", linestyle="--", linewidth=1.2,
                             label="flat-facet count")
                axis.bar(tilts, rows.gained_count, width=2.7, color="#2ca25f", alpha=0.8, label="gained vs flat")
                axis.bar(tilts, -rows.lost_count, width=2.7, color="#de2d26", alpha=0.8, label="lost vs flat")
                axis.plot(tilts, rows.total_count, "ko-", linewidth=1.4, markersize=3.5, label="total count")
                axis.set(title=f"{route} · {outcome_title} · $\\theta_0={theta0:g}^\\circ$",
                         xlabel="Fixed facet tilt [deg]", ylabel="Paired event count")
                axis.grid(alpha=0.25)
                axis.text(0.03, 0.04, "Exact streaming event-ID comparison\nrelative to the flat-facet draw.",
                          transform=axis.transAxes, fontsize=7.2, va="bottom")
        axes[0, 0].legend(fontsize=7, loc="best")
        figure.suptitle(f"Matched-event changes as fixed facet tilt increases · $\\theta_0={theta0:g}^\\circ$",
                         fontsize=13)
        figure.tight_layout()
        figure_output = pairing_output_for_boresight(output, theta0, len(boresights) > 1)
        figure_output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(figure_output, dpi=180)
        plt.close(figure)
        print(f"Wrote {figure_output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=DEFAULT_DIR / "fixed_facet_slope_step1.csv")
    parser.add_argument("--histogram-input",
                        help="compact cache from summarize_fixed_facet_slope_step1.py; avoids loading raw CSV")
    parser.add_argument("--pairing-summary",
                        help="exact paired-count cache from summarize_fixed_facet_slope_step1.py")
    parser.add_argument("--output-dir", default=DEFAULT_DIR,
                        help="directory for individual angle-distribution figures")
    parser.add_argument("--layout", choices=("individual", "grid"), default="individual",
                        help="individual: one angle per figure; grid: all angles on a compact figure")
    parser.add_argument("--grid-cols", type=int, default=3,
                        help="number of columns for --layout grid (default: 3, yielding a 2x3 layout)")
    parser.add_argument("--normalizations", default="counts,density",
                        help="comma-separated outputs: counts and/or density (default: counts,density)")
    parser.add_argument("--pairing-output", default=DEFAULT_DIR / "fixed_facet_slope_step1_pairing.png")
    parser.add_argument("--skip-pairing", action="store_true",
                        help="do not make the event-ID pairing plot (recommended with very large raw tables)")
    args = parser.parse_args()
    if args.grid_cols <= 0:
        raise ValueError("--grid-cols must be positive")
    normalizations = [item.strip().lower() for item in args.normalizations.split(",") if item.strip()]
    if not normalizations or not all(item in {"counts", "density"} for item in normalizations):
        raise ValueError("--normalizations must use counts and/or density")
    output_dir = Path(args.output_dir)
    if args.histogram_input:
        cache = pd.read_csv(args.histogram_input)
        required = {"route", "facet_tilt_deg", "lpda_theta0_deg", "selection", "field",
                    "bin_left_deg", "bin_right_deg", "count"}
        missing = required - set(cache.columns)
        if missing:
            raise ValueError(f"{args.histogram_input} lacks columns: {', '.join(sorted(missing))}")
        plot_from_histogram_cache(cache, output_dir, args.layout, args.grid_cols, normalizations)
        if args.pairing_summary:
            summary = pd.read_csv(args.pairing_summary)
            required = {"route", "lpda_theta0_deg", "selection", "facet_tilt_deg", "baseline_count",
                        "retained_count", "lost_count", "gained_count", "total_count"}
            missing = required - set(summary.columns)
            if missing:
                raise ValueError(f"{args.pairing_summary} lacks columns: {', '.join(sorted(missing))}")
            membership_from_summary(summary, Path(args.pairing_output))
        elif not args.skip_pairing:
            raise ValueError("histogram caches need --pairing-summary or --skip-pairing")
        return
    frame = pd.read_csv(args.input)
    required = {"event_id", "route", "facet_tilt_deg", "facet_aspect_deg", "lpda_theta0_deg",
                "status", "triggered", "raw_exit_angle_theta_r_deg", "theta_i_deg", "theta_emit_deg",
                "psi_deg", "payload_lat_deg", "payload_lon_deg", "payload_x_km", "payload_y_km",
                "payload_z_km", "payload_zenith_at_surface_deg", "payload_elevation_signal_deg"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{args.input} lacks columns: {', '.join(sorted(missing))}")
    for route in ("direct", "reflected"):
        if route not in set(frame.route):
            continue
        route_frame = frame[frame.route == route]
        # Valid scalar propagation is independent of LPDA boresight.  Use one
        # stored boresight row rather than duplicating identical figures.
        reference_theta0 = sorted(route_frame.lpda_theta0_deg.unique())[0]
        valid_frame = route_frame[np.isclose(route_frame.lpda_theta0_deg, reference_theta0)]
        for normalization in normalizations:
            suffix = "" if normalization == "counts" else "_density"
            if args.layout == "grid":
                distribution_grid_figure(valid_frame, output_dir / f"fixed_facet_slope_step1_{route}_valid{suffix}.png",
                                         route, "valid", args.grid_cols, normalization=normalization)
            else:
                for field, field_token, value_range, xlabel in FIELDS:
                    distribution_figure(
                        valid_frame, output_dir / f"fixed_facet_slope_step1_{route}_valid_{field_token}{suffix}.png",
                        route, "valid", field, field_token, value_range, xlabel, normalization=normalization,
                    )
        for theta0 in sorted(route_frame.lpda_theta0_deg.unique()):
            triggered_frame = route_frame[np.isclose(route_frame.lpda_theta0_deg, theta0)]
            theta0_token = f"theta0_{theta0:+g}".replace("+", "p").replace("-", "m")
            for normalization in normalizations:
                suffix = "" if normalization == "counts" else "_density"
                if args.layout == "grid":
                    distribution_grid_figure(
                        triggered_frame,
                        output_dir / f"fixed_facet_slope_step1_{route}_triggered_{theta0_token}{suffix}.png",
                        route, "triggered", args.grid_cols, theta0=theta0, normalization=normalization,
                    )
                else:
                    for field, field_token, value_range, xlabel in FIELDS:
                        distribution_figure(
                            triggered_frame,
                            output_dir / f"fixed_facet_slope_step1_{route}_triggered_{theta0_token}_{field_token}{suffix}.png",
                            route, "triggered", field, field_token, value_range, xlabel, theta0=theta0,
                            normalization=normalization,
                        )
    if not args.skip_pairing:
        membership(frame, Path(args.pairing_output))


if __name__ == "__main__":
    main()
