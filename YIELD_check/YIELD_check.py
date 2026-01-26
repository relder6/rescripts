#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import os, re, sys
import csv
import boost_histogram as bh
from tqdm import tqdm
import uproot
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import get_common_run_inputs, get_data_cuts

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, selected_target_titlename, selected_target_A, selected_target_Z = get_common_run_inputs()

# -----------------------------------------------------
# Filepaths
# -----------------------------------------------------
root_directory = "/work/hallc/c-rsidis/skimfiles/pass0"

input_settings_filepath = f"../FILTER_type/{selected_target_shortname.upper()}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_runs.dat"

# -----------------------------------------------------
# Extracting Run Information
# -----------------------------------------------------
runnums = []
prescale3_factors = []
prescale4_factors = []
beam_charges = []
tracking_effs = []
livetimes = []
hmsmomentum = []
beamcurrents = []
weights = []

with open(input_settings_filepath, "r", newline = "") as infile:
    reader = csv.DictReader(infile, delimiter = "\t")
    
    for row in reader:
        runnums.append(row["#Run#"])
        prescale3_factors.append(row["Ps3"])
        prescale4_factors.append(row["Ps4"])
        beam_charges.append(row["BCM2CutCh"])
        hmsmomentum.append(row["HMSp"])
        beamcurrents.append(row["Ibeam"])
        weights.append(row["Weight"])

if selected_target_shortname in {"c", "cu", "al", "ld2", "lh2", "dummy"}:

    branches = ["H_gtr_dp", "H_cer_npeSum", "H_cal_etottracknorm"]

    cuts = get_data_cuts()

    hist = bh.Histogram(bh.axis.Regular(100, -10, 10), storage = bh.storage.Weight())

    results = []

    for idx, (runnum, ps3, ps4, weight, charge, momentum, current) in enumerate(tqdm(zip(runnums, prescale3_factors, prescale4_factors, weights, beam_charges, hmsmomentum, beamcurrents), total = len(runnums))):
        root_filepath = f"{root_directory}/skimmed_hms_coin_replay_production_{runnum}_-1.root"
        try:
            with uproot.open(root_filepath) as rootfile:
                if "T" not in rootfile:
                    print(f"No 'T' tree in run {runnum}; skipping...")
                    continue

                tree = rootfile["T"]

                missing_branches = [b for b in branches if b not in tree.keys()]
                if missing_branches:
                    print(f"Missing branches {missing_branches} not found in run {runnum}; skipping...")
                    continue

                data = tree.arrays(branches, library = "np")

                cut = ((data["H_cer_npeSum"] > cuts["H_cer_npeSum_cut"]) &
                       (data["H_cal_etottracknorm"] > cuts["H_cal_etottracknorm_cut"]) &
                       (data["H_gtr_dp"] > cuts["H_gtr_dp_min_cut"]) & (data["H_gtr_dp"] < cuts["H_gtr_dp_max_cut"]))

                delta = "H_gtr_dp"

                if float(momentum) > 0:
                    polarity = "+"
                elif float(momentum) < 0:
                    polarity = "-"

                if selected_target_shortname in {"c", "cu", "al", "ld2", "lh2", "dummy"}:

                    run_weight = np.full_like(data[delta], float(weight) * 1000 / float(charge), dtype = float)

                    hist.reset()
                
                    hist.fill(data[delta][cut], weight=run_weight[cut])

                    delta_yield = hist.sum().value
                
                    delta_yield_err = hist.sum().variance**0.5

                    tqdm.write(f"Run {runnum}: yield = {delta_yield:.6g}, error = ±{delta_yield_err:.6g}")

                    results.append({"runnum": runnum, "polarity": polarity, "current": current, "yield": delta_yield, "yield_err": delta_yield_err})

        except FileNotFoundError:
            
            print(f"Missing file for run {runnum}, skipping...")
            
df = pd.DataFrame(results)

df["runnum"] = pd.to_numeric(df["runnum"])
df["current"] = pd.to_numeric(df["current"])
df["yield"] = pd.to_numeric(df["yield"])
df["yield_err"] = pd.to_numeric(df["yield_err"])

output_dir = "CSVs"
os.makedirs(output_dir, exist_ok=True)
output_filepath = f"{output_dir}/yield_check_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}.csv"

df.to_csv(output_filepath, index=False)
                    
print(f"Saved → {output_filepath}")
                                                     
# -----------------------------------------------------
# Plotting
# -----------------------------------------------------
elec_mask = df["polarity"] == "-"
pos_mask = df["polarity"] == "+"

output_pdf_dir = "PDFs"
os.makedirs(output_dir, exist_ok=True)
output_pdf_filepath = f"{output_pdf_dir}/yield_check_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}.pdf"

with PdfPages(output_pdf_filepath) as pdf:
    elec_df = df.loc[elec_mask]
    pos_df = df.loc[pos_mask]

    def p0_fit(y, yerr):
        w = 1.0 / (yerr**2)
        p0 = np.sum(w * y) / np.sum(w)
        p0_err = np.sqrt(1.0 / np.sum(w))
        return p0, p0_err

    elec_y = elec_df["yield"].to_numpy()
    elec_yerr = elec_df["yield_err"].to_numpy()
    elec_p0, elec_p0_err = p0_fit(elec_y, elec_yerr)

    pos_y = pos_df["yield"].to_numpy()
    pos_yerr = pos_df["yield_err"].to_numpy()
    pos_p0, pos_p0_err = p0_fit(pos_y, pos_yerr)

    tol = 0.01

    elec_hi = (1 + tol) * elec_p0
    elec_lo = (1 - tol) * elec_p0

    pos_hi = (1 + tol) * pos_p0
    pos_lo = (1 - tol) * pos_p0

    x_idx = np.arange(len(elec_df))
    x_idx_pos = np.arange(len(pos_df))

    # Yield vs Run Number page 1
    fig, (ax_top, ax_bot) = plt.subplots(2,1,figsize=(8.5,11/2), gridspec_kw={"height_ratios":[1,1], "hspace": 0.35}, sharex=False)

    fig.suptitle(f"{selected_run_type.upper()} {selected_beam_pass}Pass {selected_target_titlename} Yield vs Run Number", fontsize=12, y=0.98)
    
    ax_top.errorbar(x_idx,elec_df["yield"].to_numpy(),yerr=elec_df["yield_err"].to_numpy(),fmt="o", markersize=3, color="navy", label = "Electrons")
    ax_top.set_xticks(x_idx)
    ax_top.set_xticklabels(elec_df["runnum"].astype(str), rotation=45)
    ax_top.tick_params(axis="x", labelsize=8)
    ax_top.set_ylabel("Electron Yield", fontsize=10)
    ax_top.grid(True)
    ax_top.axhline(elec_p0, linestyle="--", color="cornflowerblue", linewidth=1.5, label="p₀ fit")
    ax_top.axhspan(elec_lo, elec_hi, color="lightskyblue", alpha=0.3, label=f"±{tol*100:.1f}% Tolerance Band")
    ax_top.legend(fontsize = 6)

    ax_bot.errorbar(x_idx_pos,pos_df["yield"].to_numpy(),yerr=pos_df["yield_err"].to_numpy(),fmt="o", markersize=3, color="red", label = "Positrons")
    ax_bot.set_xticks(x_idx_pos)
    ax_bot.set_xticklabels(pos_df["runnum"].astype(str), rotation=45)
    ax_bot.tick_params(axis="x", labelsize=8)
    ax_bot.set_ylabel("Positron Yield", fontsize=10)
    ax_bot.set_xlabel("Run Number", fontsize=10)
    ax_bot.grid(True)
    ax_bot.axhline(pos_p0, linestyle="--", color="lightcoral", alpha=0.8, linewidth=1.5, label="p₀ Fit")
    ax_bot.axhspan(pos_lo, pos_hi, color="mistyrose", label=f"±{tol*100:.1f}% Tolerance Band")
    ax_bot.legend(fontsize = 6)
    
    fig.subplots_adjust(top = 0.95, bottom = 0.13)
    pdf.savefig(fig)
    plt.close(fig)

    # Yield vs Beam Current pg 2
    fig, (ax_top, ax_bot) = plt.subplots(2,1,figsize=(8.5,11/2), gridspec_kw={"height_ratios":[1,1], "hspace": 0.15}, sharex=False)

    fig.suptitle(f"{selected_run_type.upper()} {selected_beam_pass}Pass {selected_target_titlename} Yield vs Beam Current", fontsize=12, y=0.98)

    ax_top.errorbar(df.loc[elec_mask, "current"].to_numpy(), df.loc[elec_mask, "yield"].to_numpy(), yerr = df.loc[elec_mask, "yield_err"].to_numpy(), fmt = "o", markersize = 3, color = "navy", label = "Positrons")
    ax_top.tick_params(axis="x", labelsize=8)
    ax_top.set_ylabel("Electron Yield", fontsize=10)
    ax_top.grid(True)

    ax_bot.errorbar(df.loc[pos_mask, "current"].to_numpy(), df.loc[pos_mask, "yield"].to_numpy(), yerr = df.loc[pos_mask, "yield_err"].to_numpy(), fmt = "o", markersize = 3, color = "red", label = "Data")
    ax_bot.tick_params(axis="x", labelsize=8)
    ax_bot.set_xlabel("Beam Current (µA)", fontsize=10)
    ax_bot.set_ylabel("Positron Yield", fontsize=10)
    ax_bot.grid(True)

    fig.subplots_adjust(top = 0.95, bottom = 0.10)
    pdf.savefig(fig)
    plt.close(fig)

print(f"Saved → {output_pdf_filepath}")
