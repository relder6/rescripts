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
import re

# -----------------------------------------------------
# File paths
# -----------------------------------------------------
rootfile_dir = "/work/hallc/c-rsidis/skimfiles/pass0"

# -----------------------------------------------------
# Handling user inputs to determine setting
# -----------------------------------------------------
selected_type = input("Enter desired run type (default HMSDIS): ").strip().lower()
if not selected_type:
    selected_type = f"hmsdis"
    
selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ").strip()
selected_beam_pass_to_energy_prefix = {"1": "2.",
                                       "2": "4.",
                                       "3": "6.",
                                       "4": "8.",
                                       "5": "10."}
beam_prefix = selected_beam_pass_to_energy_prefix.get(selected_beam_pass)
if not beam_prefix:
    print(f"Unknown pass: {selected_beam_pass}.  Please try again.")
    exit(1)

selected_target = input("Enter desired target (options: C, Cu, Al, LD2, LH2, Dummy): ").strip().lower()
selected_target_shortcut_to_target_variable = {"al":"al","al13":"al","aluminum":"al",
                                               "c":"c","c12":"c","carbon":"c",
                                               "cu":"cu","cu29":"cu","copper":"cu",
                                               "opt1":"optics1","optics1":"optics1",
                                               "opt2":"optics2","optics2":"optics2",
                                               "d2":"ld2","ld2":"ld2",
                                               "h2":"lh2","lh2":"lh2",
                                               "hole":"hole","chole":"hole","c-hole":"hole",
                                               "dummy":"dummy","dum":"dummy",
                                               }
selected_target_shortname = selected_target_shortcut_to_target_variable.get(selected_target)
if not selected_target_shortname:
    print(f"Unknown target: {selected_target}.  Please try again.")
    exit(1)

input_settings_filepath = f"../FILTER_type/{selected_target_shortname.upper()}/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_runs.dat"

output_dir = f"{selected_target_shortname.upper()}"

# -----------------------------------------------------
# Defining branches, using uproot to put them in data frames
# -----------------------------------------------------
runnums = []
weight = []
charge = []

with open(input_settings_filepath, "r") as infile:
    next(infile) # Skipping the header line here
    for line in infile:
        parts = line.strip().split("\t")
        runnums.append(parts[0])
        charge.append(parts[12])
        weight.append(parts[19])

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

# -----------------------------------------------------
# Histogram and csv creation
# -----------------------------------------------------
hist_data = {}
bin_edges_dict = {}

for var, bins in custom_bins.items():
    axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"], underflow=True, overflow=True)
    bin_edges_dict[var] = axis.edges
    hist_data[var] = []

for i, runnum in enumerate(runnums):
     rootfile_path = f"{rootfile_dir}/skimmed_hms_coin_replay_production_{runnum}_-1.root"
     if not os.path.exists(rootfile_path):
         print(f"WARNING: Missing {rootfile_path}, skipping...")
         continue
     df = pd.DataFrame(uproot.open(rootfile_path)["T"].arrays(branches, library = "np"))
     data_cut = (df["H_gtr_dp"].between(-8, 8) & (df["H_cer_npeSum"] > 2) & (df["H_cal_etottracknorm"] > 0.8))
     df_cut = df[data_cut].copy()

     run_weight = float(weight[i])

     for var, bins in custom_bins.items():
         axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"], underflow=True, overflow=True)
         hist = bh.Histogram(axis, storage=bh.storage.Weight())
         hist.fill(df_cut[var], weight=np.full(len(df_cut[var]), run_weight))

         counts = hist.view().value
         hist_data[var].append([runnum]+counts.tolist())

for var, rows in hist_data.items():
    bin_edges = bin_edges_dict[var]
    bin_labels = [f"bin_{i}" for i in range(len(bin_edges) - 1)]
    columns = ["runnum"] + bin_labels
    hist_df = pd.DataFrame(rows, columns=columns)
    output_filename = f"{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_histo.csv"
    output_filepath = f"{output_dir}/{output_filename}"
    hist_df.to_csv(output_filepath, index=False)
    print(f"Saved {output_filename} to {output_filepath}.")
