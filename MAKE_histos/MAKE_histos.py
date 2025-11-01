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

skimmed_rootfile_dir = "/work/hallc/c-rsidis/skimfiles/pass0"
skimmed_hms_pattern = f"skimmed_hms_coin_replay_production_{runnum}_-1.root"

rootfile_path = f"{skimmed_rootfile_dir}/{skimmed_hms_pattern}"

target = "C"

# -----------------------------------------------------
# Defining branches, using uproot to put them in data frames
# -----------------------------------------------------

branches = ["H_gtr_dp",
            "H_cal_etottracknorm",
            "H_gtr_ph",
            "H_gtr_th",
            "H_gtr_x",
            "H_gtr_y",
            "H_kin_Q2",
            "H_kin_x_bj",
            "H_kin_W",
            "H_cer_npeSum"]


df = pd.DataFrame(uproot.open(rootfile_path)["T"].arrays(branches, library="np"))
print(f"Dataframe loaded: {len(df)} rows")

data_cut = (df["H_gtr_dp"].between(-8, 8) & (df["H_cer_npeSum"] > 2) & (df["H_cal_etottracknorm"] > 0.8))

df_cut = df[data_cut].copy()
print(f"Data cuts have been applied: {data_cut}")

# -----------------------------------------------------
# Binning
# -----------------------------------------------------

custom_bins = {
"H_gtr_dp": dict(binnum = 100, min = -10, max = 10),
"H_cal_etottracknorm": dict(binnum = 100, min = 0.8, max = 1.225),
"H_gtr_ph": dict(binnum = 200, min = -0.025, max = 0.025),
"H_gtr_th": dict(binnum = 100, min = -0.075, max = 0.075),
"H_gtr_x": dict(binnum = 100, min = -0.5, max = 1.25),
"H_gtr_y": dict(binnum = 100, min = -1, max = 1.25),
"H_kin_Q2": dict(binnum = 100, min = 2.7, max = 3.9),
"H_kin_x_bj": dict(binnum = 100, min = 0.2, max = 0.3),
"H_kin_W": dict(binnum = 100, min = 3.15, max = 3.4),
"H_cer_npeSum": dict(binnum = 100, min = 0, max = 20)
}

for var, bins in custom_bins.items():
    axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"], underflow=True, overflow=True)
    hist = bh.Histogram(axis, storage=bh.storage.Weight())
    hist.fill(df_cut[var])
    bin_edges = hist.axes[0].edges
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    view = hist.view()
    counts = view.value
    errors = np.sqrt(view.variance)
    rows = []
    for i in range(len(counts)):
        rows.append((bin_centers[i], counts[i], errors[i]))
    hist_df = pd.DataFrame(rows, columns=["bin_center", "counts", "errors"])
    output_filepath = f"{target}/HMS_{runnum}_{target}_{var}_histo.csv"
    hist_df.to_csv(output_filepath, index=False)
    print(f"Histogram for {var} saved to {output_filepath}") 
