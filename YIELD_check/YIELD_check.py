#!/usr/bin/env python3

import os, re, sys, csv, uproot, mplhep
import numpy as np
import pandas as pd
import boost_histogram as bh
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from tqdm import tqdm
from matplotlib.backends.backend_pdf import PdfPages

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import parse_run_type, parse_beam_pass, parse_target, get_data_cuts

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None
arg3 = sys.argv[3] if len(sys.argv) > 3 else None

selected_run_type = parse_run_type(arg1)
selected_beam_pass, beam_prefix = parse_beam_pass(arg2)
target_abbrev, target_longname, target_shortname, target_A, target_Z = parse_target(arg3)

# -----------------------------------------------------
# Filepaths
# -----------------------------------------------------
input_settings_filepath = f"../FILTER_type/{target_abbrev.upper()}/{selected_run_type}_{selected_beam_pass}pass_{target_abbrev}_runs.csv"

histo_filepath = f"../XSEC/MAKE_csvs/{target_abbrev.upper()}/{selected_run_type}_{selected_beam_pass}pass_{target_abbrev}_H_gtr_dp_histo.csv"
err_filepath = f"../XSEC/MAKE_csvs/{target_abbrev.upper()}/{selected_run_type}_{selected_beam_pass}pass_{target_abbrev}_H_gtr_dp_err.csv"

output_png_dir = "PNGs"
os.makedirs(output_png_dir, exist_ok = True)

output_csv_dir = "CSVs"
os.makedirs(output_csv_dir, exist_ok=True)
output_csv_filepath = f"{output_csv_dir}/yield_check_{selected_run_type}_{selected_beam_pass}pass_{target_abbrev}.csv"

# -----------------------------------------------------
# Extracting Run Information
# -----------------------------------------------------
df_input = pd.read_csv(input_settings_filepath)
df_input = df_input[["runnum", "ibeam"]].copy()

df_yield = pd.read_csv(histo_filepath)
df_err = pd.read_csv(err_filepath)

runinfo_cols = ["runnum", "charge", "polarity"]
bin_cols = [c for c in df_yield.columns if c not in runinfo_cols]

# Dropping monte carlo lines here,
df_yield = df_yield[df_yield["polarity"] != "mc"].copy()
df_err = df_err[df_err["polarity"] != "mc"].copy()

# Charge normalizing by row (runnum),
df_yield[bin_cols] = df_yield[bin_cols].div(df_yield["charge"], axis=0)
df_err[bin_cols]   = df_err[bin_cols].div(df_err["charge"], axis=0)
df_yield[bin_cols] = df_yield[bin_cols].mul(1000, axis=0)
df_err[bin_cols]   = df_err[bin_cols].mul(1000, axis=0)

yield_sum = df_yield[bin_cols].sum(axis=1)
err_sum = np.sqrt((df_err[bin_cols] ** 2).sum(axis=1))

df_out = pd.DataFrame({"runnum": df_yield["runnum"],
                       "target": target_abbrev,
                       "beampass": selected_beam_pass,
                       "charge": df_yield["charge"],
                       "polarity": df_yield["polarity"],
                       "yield": yield_sum,
                       "yield_err": err_sum,})

df_out = df_out.dropna()
df_input["runnum"] = df_input["runnum"].astype(int)
df_out["runnum"]   = df_out["runnum"].astype(int)
df_out = df_out.merge(df_input, on = "runnum", how = "left")

print(df_out)

df_out.to_csv(output_csv_filepath, index = False)
                    
print(f"Saved CSV → {output_csv_filepath}")
                                                     
# -----------------------------------------------------
# Plotting
# -----------------------------------------------------
elec_mask = df_out["polarity"] == "-"
pos_mask = df_out["polarity"] == "+"

elec_df = df_out.loc[elec_mask]
pos_df = df_out.loc[pos_mask]
x_idx = np.arange(len(elec_df))
x_idx_pos = np.arange(len(pos_df))

def p0_fit(y, yerr):
    w = 1.0 / (yerr**2)
    p0 = np.sum(w * y) / np.sum(w)
    p0_err = np.sqrt(1.0 / np.sum(w))
    chi2 = np.sum(((y-p0)/yerr)**2)
    ndf = len(y)-1
    chi2_ndf = chi2/ndf if ndf > 0 else np.nan
    return p0, p0_err, chi2, chi2_ndf

def weighted_rms(y, yerr, p0):
    w = 1.0 / (yerr**2)
    variance = np.sum(w * (y - p0)**2) / np.sum(w)
    return np.sqrt(variance)

# Yield stability
elec_y = elec_df["yield"].to_numpy()
elec_yerr = elec_df["yield_err"].to_numpy()
elec_p0, elec_p0_err, elec_chi2, elec_chi2_ndf = p0_fit(elec_y, elec_yerr)
elec_sigma = weighted_rms(elec_y, elec_yerr, elec_p0)
elec_sigma_percent = elec_sigma / elec_p0 * 100
elec_hi = elec_p0 + elec_sigma
elec_lo = elec_p0 - elec_sigma

pos_y = pos_df["yield"].to_numpy()
pos_yerr = pos_df["yield_err"].to_numpy()
pos_p0, pos_p0_err, pos_chi2, pos_chi2_ndf = p0_fit(pos_y, pos_yerr)
pos_sigma = weighted_rms(pos_y, pos_yerr, pos_p0)
pos_sigma_percent = pos_sigma / pos_p0 * 100
pos_hi = pos_p0 + pos_sigma
pos_lo = pos_p0 - pos_sigma

# Yield vs current
elec_current_y = df_out.loc[elec_mask, "yield"].to_numpy()
elec_current_yerr = df_out.loc[elec_mask, "yield_err"].to_numpy()
elec_current_p0, elec_current_p0_err, elec_current_chi2, elec_current_chi2_ndf = p0_fit(elec_current_y, elec_current_yerr)
elec_current_sigma = weighted_rms(elec_current_y, elec_current_yerr, elec_current_p0)
elec_current_sigma_percent = elec_current_sigma / elec_current_p0 * 100
elec_current_hi = elec_current_p0 + elec_current_sigma
elec_current_lo = elec_current_p0 - elec_current_sigma

pos_current_y = df_out.loc[pos_mask, "yield"].to_numpy()
pos_current_yerr = df_out.loc[pos_mask, "yield_err"].to_numpy()
pos_current_p0, pos_current_p0_err, pos_current_chi2, pos_current_chi2_ndf = p0_fit(pos_current_y, pos_current_yerr)
pos_current_sigma = weighted_rms(pos_current_y, pos_current_yerr, pos_current_p0)
pos_current_sigma_percent = pos_current_sigma / pos_current_p0 * 100
pos_current_hi = pos_current_p0 + pos_current_sigma
pos_current_lo = pos_current_p0 - pos_current_sigma



# Yield vs Run Number page 1
plt.style.use(mplhep.style.ROOT)
plt.rcParams.update({"figure.titlesize": 14,
                     "axes.titlesize": 12,
                     "axes.labelsize": 10,
                     "legend.fontsize": 8,
                     "xtick.labelsize": 8,
                     "ytick.labelsize": 8})
fig, (ax_top, ax_bot) = plt.subplots(2,1,figsize=(8.5,11/2), gridspec_kw={"height_ratios":[1,1], "hspace": 0.35}, sharex=False)

fig.suptitle(f"{selected_run_type.upper()} {selected_beam_pass}Pass {target_longname} Yield vs Run Number", fontsize=12, y=0.98)
    
ax_top.axhline(elec_p0, linestyle="--", color="cornflowerblue", linewidth=1.5, label="$p_0$ fit", zorder=1)
ax_top.errorbar(x_idx,elec_df["yield"].to_numpy(),yerr=elec_df["yield_err"].to_numpy(),fmt="o", markersize=3, color="navy", label = "Electrons", zorder=2)
ax_top.set_xticks(x_idx)
ax_top.set_xticklabels(elec_df["runnum"].astype(str), rotation=45)
ax_top.set_ylabel("Electron Yield")
ax_top.grid()
ax_top.axhspan(elec_lo, elec_hi, color="lightskyblue", alpha=0.3, label=rf"±{elec_sigma_percent:.1f}% (1 $\sigma$)", zorder=0)
ax_top.text(0.02, 0.95, f"$/chi^2/ndf$ = {elec_chi2_ndf: .2f}", transform = ax_top.transAxes, verticalalignment = 'top', fontsize = 9)
ax_top.legend(loc = "upper right")

ax_bot.axhline(pos_p0, linestyle="--", color="lightcoral", alpha=0.8, linewidth=1.5, label="$p_0$ fit", zorder=1)
ax_bot.errorbar(x_idx_pos,pos_df["yield"].to_numpy(),yerr=pos_df["yield_err"].to_numpy(),fmt="o", markersize=3, color="red", label = "Positrons", zorder=2)
ax_bot.set_xticks(x_idx_pos)
ax_bot.set_xticklabels(pos_df["runnum"].astype(str), rotation=45)
ax_bot.tick_params(axis="x")
ax_bot.set_ylabel("Positron Yield")
ax_bot.set_xlabel("Run Number")
ax_bot.tick_params(axis="x")
ax_bot.grid()
ax_bot.axhspan(pos_lo, pos_hi, color="mistyrose", alpha = 0.3, label = rf"$\pm${pos_sigma_percent:.1f}% (1 $\sigma$)", zorder=0)
ax_bot.text(0.02, 0.95, f"$/chi^2/ndf$ = {pos_chi2_ndf: .2f}", transform = ax_bot.transAxes, fontsize = 9, ha = "left", va = "top")
ax_bot.legend(loc = "upper right")
    
# fig.subplots_adjust(top = 0.95, bottom = 0.13)
# pdf.savefig(fig)
fig.savefig(f"{output_png_dir}/YIELD_vs_run/yield_vs_run_{selected_run_type}_{selected_beam_pass}pass_{target_abbrev}.png", dpi=300)
plt.close(fig)

# Yield vs Beam Current pg 2
plt.style.use(mplhep.style.ROOT)
plt.rcParams.update({"figure.titlesize": 14,
                     "axes.titlesize": 12,
                     "axes.labelsize": 10,
                     "legend.fontsize": 8,
                     "xtick.labelsize": 8,
                     "ytick.labelsize": 8})
fig, (ax_top, ax_bot) = plt.subplots(2,1,figsize=(8.5,11/2), gridspec_kw={"height_ratios":[1,1], "hspace": 0.15}, sharex=False)

fig.suptitle(f"{selected_run_type.upper()} {selected_beam_pass}Pass {target_longname} Yield vs Beam Current", y=0.98)

ax_top.axhline(elec_current_p0, linestyle="--", color="cornflowerblue", linewidth=1.5, label="$p_0$ fit", zorder=1)
ax_top.errorbar(df_out.loc[elec_mask, "ibeam"].to_numpy(), df_out.loc[elec_mask, "yield"].to_numpy(), yerr = df_out.loc[elec_mask, "yield_err"].to_numpy(), fmt = "o", markersize = 3, color = "navy", label = "Positrons", zorder = 2)
ax_top.tick_params(axis="x", labelsize=8)
ax_top.axhspan(elec_current_lo, elec_current_hi, color="lightskyblue", alpha=0.3, label=rf"$\pm${elec_current_sigma_percent:.1f}% (1 $\sigma$)", zorder=0)
ax_top.set_ylabel("Electron Yield")
ax_top.text(0.02, 0.95, f"$\chi^2/ndf$ = {elec_current_chi2_ndf: .2f}", transform = ax_top.transAxes, fontsize = 9, verticalalignment = 'top')
ax_top.grid()
ax_top.legend(loc = "upper right")

ax_bot.axhline(pos_current_p0, linestyle="--", color="lightcoral", alpha=0.8, linewidth=1.5, label="$p_0$ fit", zorder=1)
ax_bot.errorbar(df_out.loc[pos_mask, "ibeam"].to_numpy(), df_out.loc[pos_mask, "yield"].to_numpy(), yerr = df_out.loc[pos_mask, "yield_err"].to_numpy(), fmt = "o", markersize = 3, color = "red", label = "Data", zorder=2)
ax_bot.tick_params(axis="x")
ax_bot.set_xlabel("Beam Current (µA)")
ax_bot.set_ylabel("Positron Yield")
ax_bot.axhspan(pos_current_lo, pos_current_hi, color="mistyrose", alpha = 0.3, label = rf"$\pm${pos_current_sigma_percent:.1f}% (1 $\sigma$)", zorder=0)
ax_bot.text(0.02, 0.95, f"$\chi^2/ndf$ = {pos_current_chi2_ndf: .2f}", transform = ax_bot.transAxes, fontsize = 9, ha = "left", va = "top")
ax_bot.grid()
ax_bot.legend(loc = "upper right")

# fig.subplots_adjust(top = 0.95, bottom = 0.10)
# pdf.savefig(fig)
fig.savefig(f"{output_png_dir}/YIELD_vs_current/yield_vs_current_{selected_run_type}_{selected_beam_pass}pass_{target_abbrev}.png", dpi=300)
plt.close(fig)

print(f"Saved PNGs → {output_png_dir}")
