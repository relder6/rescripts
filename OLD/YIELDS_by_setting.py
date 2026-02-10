#!/usr/bin/env python3

import re
import pandas as pd
import uproot
import boost_histogram as bh
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from tqdm import tqdm
import csv

# ===============================
# Filepath inputs
# ===============================

input_root_directory = "/work/hallc/c-rsidis/skimfiles/pass0"

# ===============================
# User inputs, input processing.  Will grab input file from ../FILTER_type/ directories.
# ===============================

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

# ===============================
# Extracting information
# ===============================

#Run#\tDate\ttStart\tEbeam\tIbeam\tTarget\tHMSp\tHMSth\tSHMSp\tSHMSth\tPrescaleSettings\tRunType\tBCM2CutCh\tPs3\tPs4\ttLive\tPTrigs\tELREAL\tEff\tWeight\tnu\tQ2\tepsilon\txbj\tfan_mean\tboilcorr\t# Comments

runnums = []
prescale3_factors = []
prescale4_factors = []
beam_charge = []
tracking_effs = []
livetimes = []
hmsmomentum = []
weights = []

with open(input_settings_filepath, "r", newline = "") as infile:
    reader = csv.DictReader(infile, delimiter = "\t")
    
    for row in reader:
        runnums.append(row["#Run#"])
        prescale3_factors.append(row["Ps3"])
        prescale4_factors.append(row["Ps4"])
        beam_charge.append(row["BCM2CutCh"])
        hmsmomentum.append(row["HMSp"])
        weights.append(row["Weight"])

# ===============================
# Initializing histogram information
# ===============================
        
branches = ["H_gtr_dp", "H_cer_npeSum", "H_cal_etottracknorm"]

delta = "H_gtr_dp"

hist = bh.Histogram(bh.axis.Regular(100, -10, 10), storage = bh.storage.Weight())

results = []

print(f"Starting analysis for {len(runnums)} runs...")

for idx, (runnum, ps3, ps4, charge, weight) in enumerate(tqdm(zip(runnums, prescale3_factors, prescale4_factors, beam_charge, weights), total=len(runnums))):
    ps3_val = float(ps3)
    ps4_val = float(ps4)

    if ps3_val == -999 or ps4_val == -999:
        print(f"WARNING: run {runnum} has issues with prescale values.  Fix the lookup table for this setting!")
        continue
    elif ps3_val != -1 and ps4_val != -1:
        print(f"WARNING: run {runnum} has both prescales set (ps3={ps3_val}, ps4={ps4_val}), using ps3.")
    elif ps3_val != -1 and ps4_val == -1:
        ps = ps3_val
    elif ps4_val != -1 and ps3_val == -1:
        ps = ps4_val
    else:
        tqdm.write(f"Run {runnum} has no valid prescale (both -1), skipping...")
        continue
    input_root_filepath = f"{input_root_directory}/skimmed_hms_coin_replay_production_{runnum}_-1.root"
    try:
        with uproot.open(input_root_filepath) as file:
            # print(f"Starting analysis for run {runnum}, continuing...")
            if "T" not in file:
                print(f"No 'T' tree in run {runnum}, skipping...")
                continue

            tree = file["T"]

            if delta not in tree.keys():
                print(f"Branch {delta} not found in run {runnum}, skipping...")
                continue
            
            data = tree.arrays(branches, library = "np")

            cut = (
                (data["H_cer_npeSum"] > 2) &
                (data["H_cal_etottracknorm"] > 0.8) &
                (data["H_gtr_dp"] > -8.0) &
                (data["H_gtr_dp"] < 8.0)
                )

            # Applying per run weights & charge normalizations, and converting to yield/mC
            run_weight = np.full_like(data[delta], float(weight) * 1000 / float(charge), dtype = float)

            hist.reset()
            hist.fill(data[delta][cut], weight=run_weight[cut])

            delta_yield = hist.sum().value
            delta_yield_err = hist.sum().variance**0.5

            tqdm.write(f"Run {runnum}: yield = {delta_yield:.6g}, error = ±{delta_yield_err:.6g}")

            results.append({"runnum": runnum, "yield": delta_yield, "yield_err": delta_yield_err})

            # print(f"Filled histogram for run {runnum}, moving on...")

    except FileNotFoundError:
        print(f"Missing file for run {runnum}, skipping...")
            
df = pd.DataFrame(results)

# #################

# -- Categorizing runs based on central momentum polarity for plotting with different shapes
momentum_values = np.array([float(m) for m in hmsmomentum[:len(df)]])
pos_mask = momentum_values > 0
elec_mask = momentum_values <= 0

parsed_ps3 = np.array([float(ps3) for ps3 in prescale3_factors[:len(df)]])
x_indices = np.arange(len(df))

# -- Categorizing runs based on prescale for plotting with different colors
ps3_mask = (parsed_ps3 != -999) & (parsed_ps3 != -1)
ps4_mask = ~ps3_mask

def plot_yields(df, x_indices, mask_group, title_suffix, filename_suffix):
    plt.figure(figsize=(12,6))

    # Ps3
    mask = mask_group & ps3_mask
    plt.scatter(np.array(x_indices)[mask], np.array(df.loc[mask, "yield"]),
                color='red', s=20, marker='o', label='Ps3 HMS 3/4')
    plt.errorbar(np.array(x_indices)[mask], np.array(df.loc[mask, "yield"]),
                 yerr=np.array(df.loc[mask, "yield_err"]), fmt='none', ecolor='red')

    # Ps4
    mask = mask_group & ps4_mask
    plt.scatter(np.array(x_indices)[mask], np.array(df.loc[mask, "yield"]),
                color='navy', s=20, marker='o', label='Ps4 HMS ELREAL')
    plt.errorbar(np.array(x_indices)[mask], np.array(df.loc[mask, "yield"]),
                 yerr=np.array(df.loc[mask, "yield_err"]), fmt='none', ecolor='navy')

    plt.xticks(x_indices, df["runnum"], rotation=45, fontsize=10)
    plt.xlabel("Run Number", fontsize=12)
    plt.ylabel("Delta Yield per mC", fontsize=12)
    plt.title(f"{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_{title_suffix}_yields", fontsize=14)
    plt.grid(True, linestyle='--', color='gray', alpha=0.3)
    plt.margins(y=0.1)
    plt.legend(frameon=True, fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_{filename_suffix}_yields.png")
    plt.close()

if df.empty:
    print(f"No valid runs were found for {selected_type} {selected_beam_pass}pass {selected_target_shortname}.")
    print("Verify this setting exists before trying again.  Skipping...")
    exit(0)

# Plot electrons (negative momentum)
plot_yields(df, x_indices, elec_mask, "elec", "elec")

# Plot positrons (positive momentum)
plot_yields(df, x_indices, pos_mask, "pos", "pos")
