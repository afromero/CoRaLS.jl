#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import glob
import os

# Maps user argument to dataframe column and x-label
XCOLUMNS = {
    "depth": ("Ice Depth (m)", "Ice Depth (m)"),
    "alt":   ("Altitude (km)", "Altitude (km)"),
    "ang":   ("Angle (deg)", "Angle (deg)"),
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot event rates vs ice depth or altitude, colored by energy"
    )
    parser.add_argument(
        "--input_path",
        type=str,
        default="out/",
        help="Directory where .out CSVs are stored"
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="*.out",
        help="Glob pattern for output files"
    )
    parser.add_argument(
        "--quantity",
        type=str,
        default="ARW Reflected count",
        choices=["Reflected Count", "Direct Count", "ARW Reflected count"],
        help="Which event rate to plot"
    )
    parser.add_argument(
        "--xcol",
        type=str,
        default="depth",
        choices=XCOLUMNS.keys(),
        help="Which variable to plot on the x-axis: 'depth' for ice depth, 'alt' for altitude"
    )
    parser.add_argument(
        "--title",
        type=str,
        default=None,
        help="Plot title (auto-generated if not provided)"
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show plot interactively"
    )
    parser.add_argument(
        "--save",
        type=str,
        default=None,
        help="Save plot to file (provide filename, e.g., plot.png)"
    )
    return parser.parse_args()

def get_header_info(file_path):
    with open(file_path) as fh:
        lines = [next(fh) for _ in range(6)]
        info_line = lines[1].strip()
        ntrials_line = lines[2].strip()
        return info_line, ntrials_line

def main():
    args = parse_args()
    xcol, xlabel = XCOLUMNS[args.xcol]
    
    # Gather files
    search_pattern = os.path.join(args.input_path, args.pattern)
    all_files = glob.glob(search_pattern)
    if not all_files:
        raise RuntimeError(f"No files found for pattern: {search_pattern}")

    # Read header info from the first file only
    info_line, ntrials_line = get_header_info(all_files[0])

    # Read and concatenate data from all files
    dfs = [pd.read_csv(f, header=None, skiprows=6) for f in all_files]
    df = pd.concat(dfs, ignore_index=True)
    df.columns = [
        "Energy", "Angle (deg)", "Altitude (km)", "Ice Depth (m)",
        "Reflected Count", "Reflected Error",
        "Direct Count", "Direct Error",
        "ARW Reflected count", "ARW Reflected Error"
    ]

    # Sort energies for coloring
    energies = sorted(df["Energy"].unique())
    nE = len(energies)
    cmap = mcolors.LinearSegmentedColormap.from_list("energy_cm", ["orange", "purple"], N=nE)
    colors = [cmap(i/(nE-1)) for i in range(nE)]

    fig, ax = plt.subplots(figsize=(10,6))

    # Plot per-energy lines
    '''for idx, E in enumerate(energies):
        sub = df[df["Energy"] == E].sort_values(xcol)
        ax.errorbar(
            np.array(sub[xcol]),
            2*np.array(sub[args.quantity]),
            yerr=2*np.array(sub[args.quantity.replace("Count", "Error")]),
            color=np.array(colors[idx]),
            marker='o',
            linestyle='-',
            label=f"{E:.3f}",
            alpha=0.6,
        )'''

    # Plot total
    total = df.groupby(xcol).agg({
        args.quantity: "sum",
        args.quantity.replace("count", "Error"): lambda errs: np.sqrt((errs**2).sum())
    }).reset_index()

    ax.errorbar(
        np.array(total[xcol]),
        2*np.array(total[args.quantity]),
        yerr=2*np.array(total[args.quantity.replace("count", "Error")]),
        color="black",
        marker="s",
        linestyle="--",
        linewidth=2,
        label="Total",
        alpha=0.7
    )

    print(2*np.array(total[args.quantity]))
    print(2*np.array(total[args.quantity.replace("count", "Error")]))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(args.quantity)
    if args.title is not None:
        ax.set_title(args.title)
    else:
        ax.set_title(f"{args.quantity} vs {xlabel}")

    ax.grid()
    ax.legend(title="Energy (EeV)", ncol=3, fancybox=True, loc="best")

    # Add header info as a textbox in the upper right
    info_str = ntrials_line #info_line #+ "\n" + ntrials_line
    print(info_str)
    #info_str = ntrials_line[25:-9] #info_line #+ "\n" + ntrials_line
    ax.text(
        0.64, 0.99,
        info_str,
        transform=ax.transAxes,
        fontsize=9,
        va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#f9f9f9", edgecolor="#666666", alpha=0.9),
    )

    plt.tight_layout()
    if args.save:
        plt.savefig(args.save)
    if args.show or not args.save:
        plt.show()

if __name__ == "__main__":
    main()

