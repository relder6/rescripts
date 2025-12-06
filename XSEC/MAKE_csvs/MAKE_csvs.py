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
mc_dir = "/work/hallc/c-rsidis/relder/mc-single-arm"

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

input_settings_filepath = f"../../FILTER_type/{selected_target_shortname.upper()}/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_runs.dat"

output_dir = f"{selected_target_shortname.upper()}"

# -----------------------------------------------------
# Reading in monte-carlo report to obtain normfac
# -----------------------------------------------------
if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
    mc_filepath = f"{mc_dir}/worksim/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}.root"
    if not os.path.exists(mc_filepath):
        print(f"WARNING:\tNo mc-single-arm generated root file found for {selected_type}_{selected_beam_pass}_{selected_target_shortname}. Exiting...")
        exit(1)
    mc_report_filepath = f"{mc_dir}/outfiles/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}.out"
    with open(mc_report_filepath, "r") as infile:
        normfac_line = [line for line in infile if "NORMFAC" in line.upper()]
    normfac = float(normfac_line[0].split(":")[1].split()[0]) if normfac_line else None

# -----------------------------------------------------
# Defining branches, using uproot to put them in data frames
# -----------------------------------------------------
runnums = []
weight = []
charge = []
polarity = []

with open(input_settings_filepath, "r") as infile:
    next(infile) # Skipping the header line here
    for line in infile:
        parts = line.strip().split("\t")
        runnums.append(parts[0])
        charge.append(parts[12])
        weight.append(parts[19])
        polarity.append(parts[6].strip()[0])

# -----------------------------------------------------
# Defining branches, using uproot to put them in data frames
# -----------------------------------------------------
branches = ["H_gtr_dp", "H_cal_etottracknorm", "H_gtr_ph",
            "H_gtr_th", "H_gtr_x", "H_gtr_y",
            "H_kin_Q2", "H_kin_x_bj", "H_kin_W",
            "H_cer_npeSum", "H_gtr_p"]

branches_mc = ["hsdelta", "hsyptar", "hsxptar", "q2", "xb", "w", "weight", "eprime"]

variable_mc_map = {
    "H_gtr_dp": "hsdelta",
    "H_gtr_ph": "hsyptar",
    "H_gtr_th": "hsxptar",
    "H_kin_Q2": "q2",
    "H_kin_x_bj": "xb",
    "H_kin_W": "w",
    "H_gtr_p": "eprime"
    }

# -----------------------------------------------------
# Binning
# -----------------------------------------------------
if selected_beam_pass == "4":
    custom_bins = {
        "H_gtr_dp": dict(binnum = 100, min = -10.000, max = 10.000),
        "H_gtr_ph": dict(binnum = 100, min = -0.050, max = 0.050),
        "H_gtr_th": dict(binnum = 100, min = -0.100, max = 0.100),
        "H_kin_Q2": dict(binnum = 100, min = 2.400, max = 4.200),
        "H_kin_x_bj": dict(binnum = 100, min = 0.175, max = 0.325),
        "H_kin_W": dict(binnum = 100, min = 3.100, max = 3.500),
        "H_gtr_p": dict(binnum = 151, min = 1.195, max = 2.705)
    }
    
if selected_beam_pass == "5":
    custom_bins = {
        "H_gtr_dp": dict(binnum = 100, min = -10.000, max = 10.000),
        "H_gtr_ph": dict(binnum = 100, min = -0.050, max = 0.050),
        "H_gtr_th": dict(binnum = 100, min = -0.100, max = 0.100),
        "H_kin_Q2": dict(binnum = 100, min = 2.400, max = 4.200),
        "H_kin_x_bj": dict(binnum = 100, min = 0.175, max = 0.325),
        "H_kin_W": dict(binnum = 100, min = 3.100, max = 3.500),
        "H_gtr_p": dict(binnum = 351, min = 2.995, max = 6.505)
    }

# -----------------------------------------------------
# Data histogram and csv creation
# -----------------------------------------------------
hist_data = {}
hist_err_data = {}
bin_edges_dict = {}

for var, bins in custom_bins.items():
    axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"], underflow=True, overflow=True)
    bin_edges_dict[var] = axis.edges
    hist_data[var] = []
    hist_err_data[var] = []

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
         errors = np.sqrt(hist.view().variance)
         
         hist_data[var].append([runnum, charge[i], polarity[i]]+counts.tolist())
         hist_err_data[var].append([runnum, charge[i], polarity[i]]+errors.tolist())

# -----------------------------------------------------
# Monte carlo histogram and csv creation
# -----------------------------------------------------
if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
    mc_file = uproot.open(mc_filepath)
    mc_tree = mc_file["h10"]
    df_mc = pd.DataFrame(mc_tree.arrays(branches_mc, library="np"))
    mc_cut = (df_mc["hsdelta"].between(-8, 8))
    df_mc_cut = df_mc[mc_cut].copy()
    print("DEBUG: df_mc shape:", df_mc.shape)
    print("DEBUG: df_mc columns:", df_mc.columns.tolist())
    print("DEBUG: hsdelta stats:",
          "min=", np.min(df_mc["hsdelta"]),
          "max=", np.max(df_mc["hsdelta"]))

    print("DEBUG: df_mc_cut shape:", df_mc_cut.shape)
    if df_mc_cut.empty:
        print("ERROR: df_mc_cut is EMPTY — hsdelta cut removed everything.")


    mc_hist_data = {}
    mc_hist_err = {}

for var, bins in custom_bins.items():
    mc_var = None
    
    if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        mc_var = variable_mc_map.get(var)
        
    if mc_var is None or selected_target_shortname in {"dummy", "optics1", "optics2", "hole"}:
        continue
    
    axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"])
    if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        hist_mc = bh.Histogram(axis, storage=bh.storage.Weight())

        if "weight" in df_mc_cut.columns:
            event_weights = df_mc_cut["weight"].values * normfac
        else:
            print("Weights branch not found.  Exiting...")
            print(df_mc_cut.columns.tolist())
            exit(1)
        print(f"\nDEBUG: var={var}, mc_var={mc_var}")
        print("DEBUG: data min/max:", df_mc_cut[mc_var].min(), df_mc_cut[mc_var].max())
        print("DEBUG: bins for this variable:", bins["min"], bins["max"])
        print("DEBUG: weight min/max/sum:",event_weights.min(),event_weights.max(), event_weights.sum())
        hist_mc.fill(df_mc_cut[mc_var].values, weight=event_weights)
        print("DEBUG: hist sum after fill:", hist_mc.sum())
        print("DEBUG: first 10 bin contents:", hist_mc.view().value[:10])

        mc_hist_data[var] = hist_mc.view().value.tolist()
        mc_hist_err[var] = np.sqrt(hist_mc.view().variance).tolist()


for var, rows in hist_data.items():
    bin_edges = bin_edges_dict[var]
    bin_left_edges = bin_edges[:-1]
    bin_labels = [f"{i:.6f}" for i in bin_left_edges]
    columns = ["runnum", "charge", "polarity"] + bin_labels

    # Saving here for counts
    hist_df = pd.DataFrame(rows, columns=columns)
    if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        mc_row = ["MC", "0", "mc"] + mc_hist_data[var]
        hist_df.loc[0] = mc_row
        hist_df.index = range(len(hist_df))
    else:
        hist_df.index = range(len(hist_df))
    
    output_filename = f"{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_histo.csv"
    output_filepath = f"{output_dir}/{output_filename}"
    hist_df.to_csv(output_filepath, index=False)

    # Saving here for error bars
    hist_err_df = pd.DataFrame(hist_err_data[var], columns=columns)
    if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        mc_err_row = ["MC", "0", "0"] + mc_hist_err[var]
        hist_err_df.loc[0] = mc_err_row
   
    hist_err_df.index = range(len(hist_df))
    
    output_err_filename = f"{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_err.csv"
    output_err_filepath = f"{output_dir}/{output_err_filename}"
    hist_err_df.to_csv(output_err_filepath, index=False)
    
    print(f"Saved {output_filename} and {output_err_filename} to {output_filepath}.")
