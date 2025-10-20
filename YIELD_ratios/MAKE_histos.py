#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import boost_histogram as bh
import os

runnum = input(f"Input run number: ")

# -----------------------------------------------------
# File paths
# -----------------------------------------------------

directory_filepath = "/volatile/hallc/c-rsidis/cmorean/replay_pass0a/ROOTfiles"
hms_pattern = f"hms_coin_replay_production_{runnum}_-1.root"

root_filepath = f"{directory_filepath}/{hms_pattern}"

# -----------------------------------------------------
# Defining branches, using uproot to put them in data frames
# -----------------------------------------------------

branches = ["H.gtr.dp", "H.cal.etottracknorm", "H.gtr.ph", "H.gtr.th", "H.gtr.x", "H.gtr.y",
    "H.kin.Q2", "H.kin.x_bj", "H.kin.W", "H.cer.npeSum"]

df = pd.DataFrame(uproot.open(root_filepath)["T"].arrays(branches, library="np"))
print(f"Dataframe loaded: {len(df)} rows")

data_cut = (df["H.gtr.dp"].between(-8, 8) & (df["H.cer.npeSum"] > 2) & (df["H.cal.etottracknorm"] > 0.8))

df_cut = df[data_cut].copy()
print(f"Data cuts have been applied: {data_cut}")

# -----------------------------------------------------
# Binning
# -----------------------------------------------------

custom_bins = {
"H.gtr.dp": dict(binnum = 100, min = -10, max = 10),
"H.cal.etottracknorm": dict(binnum = 100, min = 0.8, max = 1.225),
"H.gtr.ph": dict(binnum = 200, min = -0.025, max = 0.025),
"H.gtr.th": dict(binnum = 100, min = -0.075, max = 0.075),
"H.gtr.x": dict(binnum = 100, min = -0.5, max = 1.25),
"H.gtr.y": dict(binnum = 100, min = -1, max = 1.25),
"H.kin.Q2": dict(binnum = 100, min = 2.7, max = 3.9),
"H.kin.x_bj": dict(binnum = 100, min = 0.2, max = 0.3),
"H.kin.W": dict(binnum = 100, min = 3.15, max = 3.4),
"H.cer.npeSum": dict(binnum = 100, min = 0, max = 20)
}

for var, bins in custom_bins.items():
    axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"], underflow=True, overflow=True)
    hist = bh.Histogram(axis, storage=bh.storage.Weight())
    hist.fill(df_cut[var])
    bin_edges = hist.axes[0].edges
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    counts = hist.view().values()
    errors = np.sqrt(hist.view().variances())
    rows = []
    for i in range(len(counts)):
        rows.append((bin_centers[i], counts[i], errors[i]))
    hist_df = pd.DataFrame(rows, columns=["bin_center", "counts", "errors"])
    output_dir = "1d_elt2"
    os.makedirs(output_dir, exist_ok=True)
    output_filepath = os.path.join(output_dir, f"../HISTOGRAM/HMS_run_{runnum}_{var}_histo.csv")
    hist_df.to_csv(output_filepath, index=False)
    print(f"Histogram for {var} saved to {output_filepath}") 
