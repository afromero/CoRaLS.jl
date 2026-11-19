#!/usr/bin/env python3
"""Plot gained/lost adjacent-tilt distributions from the compact transition cache."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FIELDS = [
    ("theta_r", (0, 95), r"Local exit angle $\theta_r$ [deg]"),
    ("theta_i", (0, 50), r"Surface-side $\theta_i$ [deg]"),
    ("theta_emit", (0, 50), r"Ray angle at scalar source [deg]"),
    ("theta_ice", (0, 60), r"R path angle at ice $\theta_{\rm ice}$ [deg]"),
    ("psi", (0, 180), r"Shower off-axis angle $\psi$ [deg]"),
    ("payload_elevation", (-90, 0), r"Payload elevation [deg]"),
    ("payload_zenith", (0, 95), r"Payload zenith at surface [deg]"),
]


def draw(axis, panel, change, color, label, normalization):
    rows = panel[panel.change == change].sort_values("bin_left_deg")
    if rows.empty:
        return
    count = rows["count"].to_numpy(dtype=float)
    edges = np.r_[rows.bin_left_deg.to_numpy(), rows.bin_right_deg.iloc[-1]]
    if normalization == "density":
        integral = count.sum() * np.diff(edges)[0]
        if integral <= 0:
            return
        count /= integral
    axis.stairs(count, edges, color=color, linewidth=1.5, alpha=0.92, label=label)
    axis.stairs(count, edges, color=color, fill=True, alpha=0.13)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",
                        default="data/terrain_diagnostics/fixed_facet_slope_step1_transitions_histograms.csv")
    parser.add_argument("--output-dir", default="data/terrain_diagnostics/fixed_facet_transitions")
    parser.add_argument("--normalization", choices=("density", "counts"), default="density")
    parser.add_argument("--route", choices=("direct", "reflected"), action="append")
    parser.add_argument("--selection", choices=("valid", "triggered"), action="append")
    parser.add_argument("--grid-cols", type=int, default=3)
    args = parser.parse_args()
    if args.grid_cols <= 0:
        raise ValueError("--grid-cols must be positive")

    cache = pd.read_csv(args.input)
    required = {"route", "lpda_theta0_deg", "selection", "tilt_from_deg", "tilt_to_deg",
                "change", "field", "bin_left_deg", "bin_right_deg", "count"}
    missing = required - set(cache.columns)
    if missing:
        raise ValueError(f"{args.input} lacks columns: {', '.join(sorted(missing))}")
    routes = args.route or [route for route in ("direct", "reflected") if route in set(cache.route)]
    selections = args.selection or [selection for selection in ("valid", "triggered")
                                    if selection in set(cache.selection)]
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    panels = cache[["route", "lpda_theta0_deg", "selection", "tilt_from_deg", "tilt_to_deg"]].drop_duplicates()
    for row in panels.itertuples(index=False):
        if row.route not in routes or row.selection not in selections:
            continue
        panel = cache[
            (cache.route == row.route)
            & np.isclose(cache.lpda_theta0_deg, row.lpda_theta0_deg)
            & (cache.selection == row.selection)
            & np.isclose(cache.tilt_from_deg, row.tilt_from_deg)
            & np.isclose(cache.tilt_to_deg, row.tilt_to_deg)
        ]
        nrows = int(np.ceil(len(FIELDS) / args.grid_cols))
        figure, axes = plt.subplots(nrows, args.grid_cols, figsize=(7.1 * args.grid_cols, 4.35 * nrows),
                                    squeeze=False)
        for axis, (field, value_range, label) in zip(axes.ravel(), FIELDS):
            field_panel = panel[panel.field == field]
            draw(axis, field_panel, "lost", "#de2d26",
                 fr"lost: evaluated at {row.tilt_from_deg:g}°", args.normalization)
            draw(axis, field_panel, "gained", "#2ca25f",
                 fr"gained: evaluated at {row.tilt_to_deg:g}°", args.normalization)
            axis.set(title=label, xlabel=label,
                     ylabel=("Normalized density [deg$^{-1}$]" if args.normalization == "density"
                             else "Event count per angle bin"),
                     xlim=value_range)
            axis.grid(alpha=0.25)
            if field == "theta_r":
                axis.axvline(90, color="0.25", linestyle=":", linewidth=1)
        for axis in axes.ravel()[len(FIELDS):]:
            axis.set_visible(False)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            axes[0, 0].legend(handles, labels, fontsize=8, loc="best")
        figure.suptitle(
            f"{row.route} · {row.selection} · facet tilt {row.tilt_from_deg:g}° → {row.tilt_to_deg:g}°"
            + "\nLost events use their last passing state; gained events use their first passing state.",
            fontsize=13,
        )
        figure.tight_layout()
        theta_token = f"theta0_{row.lpda_theta0_deg:+g}".replace("+", "p").replace("-", "m")
        filename = (f"fixed_facet_transition_{row.route}_{row.selection}_{theta_token}_"
                    f"{row.tilt_from_deg:g}_to_{row.tilt_to_deg:g}_{args.normalization}.png")
        output = output_dir / filename
        figure.savefig(output, dpi=180)
        plt.close(figure)
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
