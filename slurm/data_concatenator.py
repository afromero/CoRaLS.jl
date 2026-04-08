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
import datetime

## arguments
parser = argparse.ArgumentParser()
parser.add_argument(
        "--input_files", 
        nargs = "+",
        type = str, 
        default = ["out/*.out"], 
        help = "Directory containing .out files from a CoRaLS run."
        )
parser.add_argument(
        "--output_file",
        type = str,
        default = r'run_data_{}.csv'.format(datetime.datetime.now().timestamp()),
        help = "Output concatenated data file."
        )
args = parser.parse_args()

## Read in the csv files
dfs = [pd.read_csv(f, header = None, skiprows = 6) for f in args.input_files] 
df = pd.concat(dfs, ignore_index = True)
df.columns = [
        "Energy (EeV)", "Altitude (km)", "Ice Depth (m)", "Angle (deg)",
        "Reflected Count", "Reflected Error",
        "Direct Count", "Direct Error"
    ]
df.to_csv(args.output_file, index = False)
