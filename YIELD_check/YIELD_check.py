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
root_directory = "/w/hallc-scshelf2102/c-rsidis/skimfiles/pass0p1"

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
        chi2 = np.sum(((y-p0)/yerr)**2)
        ndof = len(y)-1
        chi2_ndof = chi2/ndof if ndof > 0 else np.nan
        return p0, p0_err, chi2, chi2_ndof

    def weighted_rms(y, yerr, p0):
        w = 1.0 / (yerr**2)
        variance = np.sum(w * (y - p0)**2) / np.sum(w)
        return np.sqrt(variance)

    #Page 1 yield stability fit, sigma
    elec_y = elec_df["yield"].to_numpy()
    elec_yerr = elec_df["yield_err"].to_numpy()
    elec_p0, elec_p0_err, elec_chi2, elec_chi2_ndof = p0_fit(elec_y, elec_yerr)
    elec_sigma = weighted_rms(elec_y, elec_yerr, elec_p0)
    elec_sigma_percent = elec_sigma / elec_p0 * 100
    elec_hi = elec_p0 + elec_sigma
    elec_lo = elec_p0 - elec_sigma

    pos_y = pos_df["yield"].to_numpy()
    pos_yerr = pos_df["yield_err"].to_numpy()
    pos_p0, pos_p0_err, pos_chi2, pos_chi2_ndof = p0_fit(pos_y, pos_yerr)
    pos_sigma = weighted_rms(pos_y, pos_yerr, pos_p0)
    pos_sigma_percent = pos_sigma / pos_p0 * 100
    pos_hi = pos_p0 + pos_sigma
    pos_lo = pos_p0 - pos_sigma

    #Page 2 yield vs current fit, sigma
    elec_current_y = df.loc[elec_mask, "yield"].to_numpy()
    elec_current_yerr = df.loc[elec_mask, "yield_err"].to_numpy()
    elec_current_p0, elec_current_p0_err, elec_current_chi2, elec_current_chi2_ndof = p0_fit(elec_current_y, elec_current_yerr)
    elec_current_sigma = weighted_rms(elec_current_y, elec_current_yerr, elec_current_p0)
    elec_current_sigma_percent = elec_current_sigma / elec_current_p0 * 100
    elec_current_hi = elec_current_p0 + elec_current_sigma
    elec_current_lo = elec_current_p0 - elec_current_sigma

    pos_current_y = df.loc[pos_mask, "yield"].to_numpy()
    pos_current_yerr = df.loc[pos_mask, "yield_err"].to_numpy()
    pos_current_p0, pos_current_p0_err, pos_current_chi2, pos_current_chi2_ndof = p0_fit(pos_current_y, pos_current_yerr)
    pos_current_sigma = weighted_rms(pos_current_y, pos_current_yerr, pos_current_p0)
    pos_current_sigma_percent = pos_current_sigma / pos_current_p0 * 100
    pos_current_hi = pos_current_p0 + pos_current_sigma
    pos_current_lo = pos_current_p0 - pos_current_sigma


    x_idx = np.arange(len(elec_df))
    x_idx_pos = np.arange(len(pos_df))

    # Yield vs Run Number page 1
    fig, (ax_top, ax_bot) = plt.subplots(2,1,figsize=(8.5,11/2), gridspec_kw={"height_ratios":[1,1], "hspace": 0.35}, sharex=False)

    fig.suptitle(f"{selected_run_type.upper()} {selected_beam_pass}Pass {selected_target_titlename} Yield vs Run Number", fontsize=12, y=0.98)
    
    ax_top.axhline(elec_p0, linestyle="--", color="cornflowerblue", linewidth=1.5, label="p₀ fit", zorder=1)
    ax_top.errorbar(x_idx,elec_df["yield"].to_numpy(),yerr=elec_df["yield_err"].to_numpy(),fmt="o", markersize=3, color="navy", label = "Electrons", zorder=2)
    ax_top.set_xticks(x_idx)
    ax_top.set_xticklabels(elec_df["runnum"].astype(str), rotation=45)
    ax_top.tick_params(axis="x", labelsize=8)
    ax_top.set_ylabel("Electron Yield", fontsize=10)
    ax_top.grid(True)
    ax_top.axhspan(elec_lo, elec_hi, color="lightskyblue", alpha=0.3, label=rf"±{elec_sigma_percent:.1f}% (1 $\sigma$)", zorder=0)
    ax_top.text(0.02, 0.95, f"$\\chi^2/ndof = {elec_chi2_ndof: .2f}$", transform = ax_top.transAxes, fontsize = 9, verticalalignment = 'top')
    ax_top.legend(fontsize = 6)

    ax_bot.axhline(pos_p0, linestyle="--", color="lightcoral", alpha=0.8, linewidth=1.5, label="p₀ Fit", zorder=1)
    ax_bot.errorbar(x_idx_pos,pos_df["yield"].to_numpy(),yerr=pos_df["yield_err"].to_numpy(),fmt="o", markersize=3, color="red", label = "Positrons", zorder=2)
    ax_bot.set_xticks(x_idx_pos)
    ax_bot.set_xticklabels(pos_df["runnum"].astype(str), rotation=45)
    ax_bot.tick_params(axis="x", labelsize=8)
    ax_bot.set_ylabel("Positron Yield", fontsize=10)
    ax_bot.set_xlabel("Run Number", fontsize=10)
    ax_bot.grid(True)
    ax_bot.axhspan(pos_lo, pos_hi, color="mistyrose", alpha = 0.3, label = rf"±{pos_sigma_percent:.1f}% (1 $\sigma$)", zorder=0)
    ax_bot.text(0.02, 0.95, f"$\\chi^2/ndof = {pos_chi2_ndof: .2f}$", transform = ax_bot.transAxes, fontsize = 9, verticalalignment = 'top')
    ax_bot.legend(fontsize = 6)
    
    fig.subplots_adjust(top = 0.95, bottom = 0.13)
    pdf.savefig(fig)
    plt.close(fig)

    # Yield vs Beam Current pg 2
    fig, (ax_top, ax_bot) = plt.subplots(2,1,figsize=(8.5,11/2), gridspec_kw={"height_ratios":[1,1], "hspace": 0.15}, sharex=False)

    fig.suptitle(f"{selected_run_type.upper()} {selected_beam_pass}Pass {selected_target_titlename} Yield vs Beam Current", fontsize=12, y=0.98)

    ax_top.axhline(elec_current_p0, linestyle="--", color="cornflowerblue", linewidth=1.5, label="p₀ fit", zorder=1)
    ax_top.errorbar(df.loc[elec_mask, "current"].to_numpy(), df.loc[elec_mask, "yield"].to_numpy(), yerr = df.loc[elec_mask, "yield_err"].to_numpy(), fmt = "o", markersize = 3, color = "navy", label = "Positrons", zorder = 2)
    ax_top.tick_params(axis="x", labelsize=8)
    ax_top.axhspan(elec_current_lo, elec_current_hi, color="lightskyblue", alpha=0.3, label=rf"±{elec_current_sigma_percent:.1f}% (1 $\sigma$)", zorder=0)
    ax_top.set_ylabel("Electron Yield", fontsize=10)
    ax_top.text(0.02, 0.95, f"$\\chi^2/ndof = {elec_current_chi2_ndof: .2f}$", transform = ax_top.transAxes, fontsize = 9, verticalalignment = 'top')
    ax_top.grid(True)

    ax_bot.axhline(pos_current_p0, linestyle="--", color="lightcoral", alpha=0.8, linewidth=1.5, label="p₀ Fit", zorder=1)
    ax_bot.errorbar(df.loc[pos_mask, "current"].to_numpy(), df.loc[pos_mask, "yield"].to_numpy(), yerr = df.loc[pos_mask, "yield_err"].to_numpy(), fmt = "o", markersize = 3, color = "red", label = "Data", zorder=2)
    ax_bot.tick_params(axis="x", labelsize=8)
    ax_bot.set_xlabel("Beam Current (µA)", fontsize=10)
    ax_bot.set_ylabel("Positron Yield", fontsize=10)
    ax_bot.axhspan(pos_current_lo, pos_current_hi, color="mistyrose", alpha = 0.3, label = rf"±{pos_current_sigma_percent:.1f}% (1 $\sigma$)", zorder=0)
    ax_bot.text(0.02, 0.95, f"$\\chi^2/ndof = {pos_current_chi2_ndof: .2f}$", transform = ax_bot.transAxes, fontsize = 9, verticalalignment = 'top')
    ax_bot.grid(True)

    fig.subplots_adjust(top = 0.95, bottom = 0.10)
    pdf.savefig(fig)
    plt.close(fig)

print(f"Saved → {output_pdf_filepath}")
