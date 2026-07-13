## File name: data_concatenator.py
## author: Alex Machtay (machtay.1@osu.edu)
## Purpose:
##      Reads outputs from a job series and concatenates them into a single file
##      This makes the data easy to contain and read in with a streamlined plotter

## Imports
import numpy as np
import pandas as pd
import glob
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# Maps user argument to dataframe column to use for plotting
COLUMNS = {
    "depth": ("Ice Depth (m)", "Ice Depth (m)"),
    "alt":   ("Altitude (km)", "Altitude (km)"),
    "ang":   ("Angle (deg)", "Angle (deg)"),
    "energy": ("Energy (EeV)", "Energy (EeV)"),
    "ref_count": ("Reflected Count", "Reflected Count"),
    "dir_count": ("Direct Count", "Direct Count"),
    "ref_error": ("Reflected Error", "Reflected Error"),
    "dir_error": ("Direct Error", "Direct Error"),
}

## arguments
parser = argparse.ArgumentParser()
parser.add_argument(
        "--input_file", 
        type = str, 
        default = "run_data.csv", 
        help = "Concatenated datafile from a series of CoRaLS runs."
        )
parser.add_argument(
        "--y_axis",
        type = str,
        default = "ref_count",
        choices = COLUMNS.keys(),
        help = "Which rate to use (direct or reflected)"
        )
parser.add_argument(
        "--x_axis",
        type = str,
        default = "alt",
        choices = COLUMNS.keys(),
        help = "Which variable for the x-axis; depth, alt, ang"
        )
parser.add_argument(
        "--line_sorting",
        type = str,
        default = "energy",
        choices = COLUMNS.keys()
        )
args = parser.parse_args()
## Read in the datafile
df = pd.read_csv(args.input_file)
line_col, _ = COLUMNS[args.line_sorting]
xcol, xlabel = COLUMNS[args.x_axis]
ycol, ylabel = COLUMNS[args.y_axis]

## We need to set the bins we are making lines of
bins = sorted(df[line_col].unique()) ## Ex: use energy
Nbins = len(bins) ## Get the length for making a color map
cmap = mcolors.LinearSegmentedColormap.from_list("energy_cm", ["orange", "purple"], N = Nbins)
colors = [cmap(i/(Nbins-1)) for i in range(Nbins)]

## Let's also include a total
x_vals = (df.groupby(xcol)[xcol].unique() )
totals = (df.groupby(xcol)[ycol].sum())
total_errs = (df.groupby(xcol)[ycol.replace("Count", "Error")].sum() )

print(len(totals))
print(len(total_errs))
print(x_vals)
print(totals)
print(total_errs)

## Make a plot
fig, ax = plt.subplots(figsize = (12, 8))
#yerr = [2*np.array(sub[args.y_axis.replace("Count", "Error")]) if COLUMNS[args.y_axis] == ("ref_count" or "dir_count") else 0]

## We make multiple lines per bin
for idx, E in enumerate(bins):
    sub = df[df[line_col] == E].sort_values(xcol)
    
    yerr_col = ycol.replace("Count", "Error")
    yerr = 2*np.array(sub[yerr_col]) if yerr_col in sub.columns else None

    ax.errorbar(
        np.array(sub[xcol]),
        2*np.array(sub[ycol]),
        yerr = yerr,
        color=np.array(colors[idx]),
        marker='o',
        linestyle='-',
        label=f"{E:.3f}",
        alpha=0.6,
    )
## Plot the total
print("Just before plotting total")
ax.errorbar(
    x_vals,
    2*totals,
    yerr = 2*total_errs,
    color='k',
    marker='o',
    linestyle='-',
    label=f"Total",
    alpha=0.6,
)



ax.legend(fontsize = 12, ncol = 3, loc='best')
#ax.set_xticks(fontsize = 14)
#ax.set_yticks(fontsize = 14)
ax.set_xlabel(xlabel, fontsize = 16)
ax.set_ylabel(ylabel, fontsize = 16)
ax.tick_params(axis = 'both', labelsize=14, which = 'major')
fig.tight_layout()
fig.savefig("new_test.png")
