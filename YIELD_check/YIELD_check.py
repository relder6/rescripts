#!/usr/bin/env python3

import os, sys, csv, uproot, mplhep
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter
from matplotlib.backends.backend_pdf import PdfPages

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)
from INIT.config import parse_run_type, parse_beam_pass, parse_target, get_data_cuts, parse_phase
plt.style.use(mplhep.style.ROOT)
plt.style.use(os.path.join(BASE_DIR, "INIT", "myhep.mplstyle"))

arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None
arg3 = sys.argv[3] if len(sys.argv) > 3 else None
arg4 = sys.argv[4] if len(sys.argv) > 4 else None

selected_run_type = parse_run_type(arg1)
selected_beam_pass, beam_prefix = parse_beam_pass(arg2)
target_abbrev, target_longname, target_shortname, target_A, target_Z = parse_target(arg3)
phase = parse_phase(arg4)

input_settings_filepath = f"../FILTER_type/{target_abbrev.upper()}/filtered_{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}.csv"
histo_filepath = f"../XSEC/MAKE_csvs/{target_abbrev.upper()}/{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}_H_gtr_dp_histo.csv"
err_filepath = f"../XSEC/MAKE_csvs/{target_abbrev.upper()}/{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}_H_gtr_dp_err.csv"

required_files = [input_settings_filepath, histo_filepath, err_filepath]
for f in required_files:
    if not os.path.exists(f):
        print(f"Required file missing: {f}")
        sys.exit(1)

output_pdf_dir = "PDFs"
os.makedirs(output_pdf_dir, exist_ok=True)

output_csv_dir = "CSVs"
os.makedirs(output_csv_dir, exist_ok=True)

output_csv_filepath = f"{output_csv_dir}/yield_check_{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}.csv"
output_scaler_csv_filepath = f"{output_csv_dir}/scaler_yields_{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}.csv"

pdf_filepath = f"{output_pdf_dir}/YIELD_check_{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}.pdf"

df_yield = pd.read_csv(histo_filepath)
df_err = pd.read_csv(err_filepath)
df_yield = df_yield[df_yield["polarity"] != "mc"].copy()
df_err = df_err[df_err["polarity"] != "mc"].copy()

df_yield.iloc[:, 3:] = df_yield.iloc[:, 3:].div(df_yield["charge"], axis=0).mul(1000)
df_err.iloc[:, 3:] = df_err.iloc[:, 3:].div(df_err["charge"], axis=0).mul(1000)

df = df_yield[["runnum", "charge", "polarity"]].copy()
df["yield"] = df_yield.iloc[:, 3:].sum(axis=1)
df["yield_err"] = np.sqrt((df_err.iloc[:, 3:] ** 2).sum(axis=1))

df_input = pd.read_csv(input_settings_filepath, usecols=["runnum", "ibeam", "qbeam", "elclean_counts", "hms_p"])
df_input = df_input.assign(
    runnum=pd.to_numeric(df_input["runnum"], errors="coerce"),
    qbeam=pd.to_numeric(df_input["qbeam"], errors="coerce"),
    ibeam=pd.to_numeric(df_input["ibeam"], errors="coerce"),
    elclean_counts=pd.to_numeric(df_input["elclean_counts"], errors="coerce"),
    hms_p=pd.to_numeric(df_input["hms_p"], errors="coerce"),
)
df_input = df_input.dropna(subset=["runnum", "qbeam", "ibeam", "elclean_counts", "hms_p"]).copy()
df_input["runnum"] = df_input["runnum"].astype(int)

df["runnum"] = pd.to_numeric(df["runnum"], errors="coerce")
df = df.dropna(subset=["runnum"]).copy()
df["runnum"] = df["runnum"].astype(int)
df = df.merge(df_input[["runnum", "ibeam"]], on="runnum", how="left")

df["target"] = target_abbrev
df["beampass"] = selected_beam_pass

df.to_csv(output_csv_filepath, index=False)
print(f"Saved tracked CSV → {output_csv_filepath}")

df_scaler = df_input.copy()
df_scaler["target"] = target_abbrev
df_scaler["beampass"] = selected_beam_pass
df_scaler["polarity"] = np.where(df_scaler["hms_p"] < 0, "-", "+")
df_scaler["scaler_yield"] = df_scaler["elclean_counts"] / df_scaler["qbeam"]
df_scaler["scaler_err"] = np.sqrt(df_scaler["elclean_counts"]) / df_scaler["qbeam"]
df_scaler.to_csv(output_scaler_csv_filepath, index=False)
print(f"Saved scalers CSV → {output_scaler_csv_filepath}")

elec_mask = df["polarity"] == "-"
pos_mask = df["polarity"] == "+"
elec_df = df.loc[elec_mask].copy()
pos_df = df.loc[pos_mask].copy()
elec_df = elec_df[np.isfinite(elec_df["yield_err"]) & (elec_df["yield_err"] > 0)].copy()
pos_df = pos_df[np.isfinite(pos_df["yield_err"]) & (pos_df["yield_err"] > 0)].copy()

x_idx = np.arange(len(elec_df))
x_idx_pos = np.arange(len(pos_df))

def p0_fit(y, yerr):
    mask = np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    y = y[mask]
    yerr = yerr[mask]
    if len(y) == 0:
        return np.nan, np.nan, np.nan, np.nan
    w = 1.0 / (yerr**2)
    p0 = np.sum(w * y) / np.sum(w)
    p0_err = np.sqrt(1.0 / np.sum(w))
    chi2 = np.sum(((y - p0) / yerr) ** 2)
    ndf = len(y) - 1
    chi2_ndf = chi2 / ndf if ndf > 0 else np.nan
    return p0, p0_err, chi2, chi2_ndf

def weighted_rms(y, yerr, p0):
    mask = np.isfinite(y) & np.isfinite(yerr) & (yerr > 0) & np.isfinite(p0)
    y = y[mask]
    yerr = yerr[mask]
    if len(y) == 0:
        return np.nan
    w = 1.0 / (yerr**2)
    variance = np.sum(w * (y - p0) ** 2) / np.sum(w)
    return np.sqrt(variance)

elec_y = elec_df["yield"].to_numpy()
elec_yerr = elec_df["yield_err"].to_numpy()
elec_p0, elec_p0_err, elec_chi2, elec_chi2_ndf = p0_fit(elec_y, elec_yerr)
elec_sigma = weighted_rms(elec_y, elec_yerr, elec_p0)
elec_sigma_percent = elec_sigma / elec_p0 * 100 if np.isfinite(elec_p0) and elec_p0 != 0 else np.nan
elec_hi = elec_p0 + elec_sigma
elec_lo = elec_p0 - elec_sigma

pos_y = pos_df["yield"].to_numpy()
pos_yerr = pos_df["yield_err"].to_numpy()
pos_p0, pos_p0_err, pos_chi2, pos_chi2_ndf = p0_fit(pos_y, pos_yerr)
pos_sigma = weighted_rms(pos_y, pos_yerr, pos_p0)
pos_sigma_percent = pos_sigma / pos_p0 * 100 if np.isfinite(pos_p0) and pos_p0 != 0 else np.nan
pos_hi = pos_p0 + pos_sigma
pos_lo = pos_p0 - pos_sigma

elec_current_y = elec_df["yield"].to_numpy()
elec_current_yerr = elec_df["yield_err"].to_numpy()
elec_current_p0, elec_current_p0_err, elec_current_chi2, elec_current_chi2_ndf = p0_fit(elec_current_y, elec_current_yerr)
elec_current_sigma = weighted_rms(elec_current_y, elec_current_yerr, elec_current_p0)
elec_current_sigma_percent = elec_current_sigma / elec_current_p0 * 100 if np.isfinite(elec_current_p0) and elec_current_p0 != 0 else np.nan
elec_current_hi = elec_current_p0 + elec_current_sigma
elec_current_lo = elec_current_p0 - elec_current_sigma

pos_current_y = pos_df["yield"].to_numpy()
pos_current_yerr = pos_df["yield_err"].to_numpy()
pos_current_p0, pos_current_p0_err, pos_current_chi2, pos_current_chi2_ndf = p0_fit(pos_current_y, pos_current_yerr)
pos_current_sigma = weighted_rms(pos_current_y, pos_current_yerr, pos_current_p0)
pos_current_sigma_percent = pos_current_sigma / pos_current_p0 * 100 if np.isfinite(pos_current_p0) and pos_current_p0 != 0 else np.nan
pos_current_hi = pos_current_p0 + pos_current_sigma
pos_current_lo = pos_current_p0 - pos_current_sigma

scaler_df = df_scaler[np.isfinite(df_scaler["scaler_yield"]) & np.isfinite(df_scaler["scaler_err"]) & (df_scaler["scaler_err"] > 0)].copy()
scaler_elec_mask = scaler_df["polarity"] == "-"
scaler_pos_mask = scaler_df["polarity"] == "+"

scaler_elec_df = scaler_df.loc[scaler_elec_mask].copy()
scaler_pos_df = scaler_df.loc[scaler_pos_mask].copy()

scaler_x_idx = np.arange(len(scaler_elec_df))
scaler_x_idx_pos = np.arange(len(scaler_pos_df))

scaler_elec_y = scaler_elec_df["scaler_yield"].to_numpy()
scaler_elec_yerr = scaler_elec_df["scaler_err"].to_numpy()
scaler_elec_p0, scaler_elec_p0_err, scaler_elec_chi2, scaler_elec_chi2_ndf = p0_fit(scaler_elec_y, scaler_elec_yerr)
scaler_elec_sigma = weighted_rms(scaler_elec_y, scaler_elec_yerr, scaler_elec_p0)
scaler_elec_sigma_percent = scaler_elec_sigma / scaler_elec_p0 * 100 if np.isfinite(scaler_elec_p0) and scaler_elec_p0 != 0 else np.nan

scaler_pos_y = scaler_pos_df["scaler_yield"].to_numpy()
scaler_pos_yerr = scaler_pos_df["scaler_err"].to_numpy()
scaler_pos_p0, scaler_pos_p0_err, scaler_pos_chi2, scaler_pos_chi2_ndf = p0_fit(scaler_pos_y, scaler_pos_yerr)
scaler_pos_sigma = weighted_rms(scaler_pos_y, scaler_pos_yerr, scaler_pos_p0)
scaler_pos_sigma_percent = scaler_pos_sigma / scaler_pos_p0 * 100 if np.isfinite(scaler_pos_p0) and scaler_pos_p0 != 0 else np.nan

scaler_elec_current_y = scaler_elec_df["scaler_yield"].to_numpy()
scaler_elec_current_yerr = scaler_elec_df["scaler_err"].to_numpy()
scaler_elec_current_p0, scaler_elec_current_p0_err, scaler_elec_current_chi2, scaler_elec_current_chi2_ndf = p0_fit(scaler_elec_current_y, scaler_elec_current_yerr)
scaler_elec_current_sigma = weighted_rms(scaler_elec_current_y, scaler_elec_current_yerr, scaler_elec_current_p0)
scaler_elec_current_sigma_percent = scaler_elec_current_sigma / scaler_elec_current_p0 * 100 if np.isfinite(scaler_elec_current_p0) and scaler_elec_current_p0 != 0 else np.nan
scaler_elec_current_hi = scaler_elec_current_p0 + scaler_elec_current_sigma
scaler_elec_current_lo = scaler_elec_current_p0 - scaler_elec_current_sigma

scaler_pos_current_y = scaler_pos_df["scaler_yield"].to_numpy()
scaler_pos_current_yerr = scaler_pos_df["scaler_err"].to_numpy()
scaler_pos_current_p0, scaler_pos_current_p0_err, scaler_pos_current_chi2, scaler_pos_current_chi2_ndf = p0_fit(scaler_pos_current_y, scaler_pos_current_yerr)
scaler_pos_current_sigma = weighted_rms(scaler_pos_current_y, scaler_pos_current_yerr, scaler_pos_current_p0)
scaler_pos_current_sigma_percent = scaler_pos_current_sigma / scaler_pos_current_p0 * 100 if np.isfinite(scaler_pos_current_p0) and scaler_pos_current_p0 != 0 else np.nan
scaler_pos_current_hi = scaler_pos_current_p0 + scaler_pos_current_sigma
scaler_pos_current_lo = scaler_pos_current_p0 - scaler_pos_current_sigma

have_positrons = not pos_df.empty and not scaler_pos_df.empty

# -----------------------------------------------------
# Plotting
# -----------------------------------------------------
# plt.style.use(mplhep.style.ROOT)
# plt.rcParams.update({"figure.figsize": (8.5, 6.5),
#                      "figure.titlesize": 20,
#                      "font.size": 16,
#                      "axes.grid": True,
#                      "grid.alpha": 0.5,
#                      "legend.frameon": True,
#                      "lines.markersize": 5,
#                      "xtick.top": False, "ytick.right": False,
#                      "xtick.direction": "in", "ytick.direction": "in",
#                      "xtick.minor.visible": True, "ytick.minor.visible": True,
#                      "xtick.major.size": 10, "ytick.major.size": 10,
#                      "xtick.major.width": 1.5, "ytick.major.width": 1.5,
#                      "xtick.minor.size": 5, "ytick.minor.size": 5,
#                      "xtick.minor.width": 1.0, "ytick.minor.width": 1.0,})

with PdfPages(pdf_filepath) as pdf:

    # -------------------------------------------------
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 1], "hspace": 0},sharex=True)

    fig.suptitle(f"{selected_run_type.upper()} {selected_beam_pass}Pass {target_longname} Electron Yield vs Run Number")

    ax_top.axhline(elec_p0, linestyle="-", color="cornflowerblue", linewidth=2.0, label="Tracked $p_0$ fit", zorder=1)
    ax_top.axhline(elec_p0*1.01, linestyle="--", color="cornflowerblue", linewidth=1.5, label="Tracked $p_0$ + 1%")
    ax_top.axhline(elec_p0*0.99, linestyle="--", color="cornflowerblue", linewidth=1.5, label="Tracked $p_0$ - 1%")
    ax_top.errorbar(x_idx, elec_df["yield"].to_numpy(), yerr=elec_df["yield_err"].to_numpy(), fmt="o",  color="navy", label="Tracked", zorder=2)
    ax_top.set_ylabel("Electron Tracked Yield")
    ax_top.text(0.02, 0.95, f"$\\chi^2/ndf$ = {elec_chi2_ndf: .2f}", transform=ax_top.transAxes, verticalalignment="top")
    ax_top.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    ax_top.tick_params(labelbottom=False)
    
    ax_bot.axhline(scaler_elec_p0, linestyle="-", color="darkred", linewidth=2.0, label="Scaler $p_0$ fit", zorder=1)
    ax_bot.axhline(scaler_elec_p0*1.01, linestyle="--", color="darkred", linewidth=1.5, label="Scaler $p_0$ + 1%")
    ax_bot.axhline(scaler_elec_p0*0.99, linestyle="--", color="darkred", linewidth=1.5, label="Scaler $p_0$ - 1%")
    ax_bot.errorbar(scaler_x_idx, scaler_elec_df["scaler_yield"].to_numpy(), yerr=scaler_elec_df["scaler_err"].to_numpy(), fmt="o", color="red", label="Scaler", zorder=2)
    ax_bot.set_xticks(x_idx)
    ax_bot.set_xticklabels(elec_df["runnum"].astype(str), rotation=45)
    ax_bot.set_ylabel("Electron Scaler Yield")
    ax_bot.set_xlabel("Run Number")
    ax_bot.text(0.02, 0.95, f"$\\chi^2/ndf$ = {scaler_elec_chi2_ndf: .2f}", transform=ax_bot.transAxes, ha="left", va="top")
    ax_bot.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    # set_min_ylim(ax_bot, scaler_elec_p0, scaler_elec_df["scaler_yield"])

    #Settings for tick marks
    ax_top.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_top.yaxis.set_minor_locator(AutoMinorLocator(5))

    ax_bot.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_bot.yaxis.set_minor_locator(AutoMinorLocator(5))
    
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    # -------------------------------------------------
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 1], "hspace": 0},sharex=True)

    fig.suptitle(f"{selected_run_type.upper()} {selected_beam_pass}Pass {target_longname} Electron Yield vs Beam Current", y=0.98)

    ax_top.axhline(elec_current_p0, linestyle="-", color="cornflowerblue", linewidth=2.0, label="Tracked $p_0$ fit", zorder=1)
    ax_top.axhline(elec_current_p0*1.01, linestyle="--", color="cornflowerblue", linewidth=1.5, label="Tracked $p_0$ + 1%")
    ax_top.axhline(elec_current_p0*0.99, linestyle="--", color="cornflowerblue", linewidth=1.5, label="Tracked $p_0$ - 1%")

    ax_top.errorbar(elec_df["ibeam"].to_numpy(), elec_df["yield"].to_numpy(), yerr=elec_df["yield_err"].to_numpy(), fmt="o",  color="navy", label="Tracked", zorder=2)
    ax_top.set_ylabel("Electron Tracked Yield")
    ax_top.text(0.02, 0.95, f"$\\chi^2/ndf$ = {elec_current_chi2_ndf: .2f}", transform=ax_top.transAxes, verticalalignment="top")
    ax_top.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    ax_top.tick_params(labelbottom=False)

    ax_bot.axhline(scaler_elec_current_p0, linestyle="-", color="darkred", linewidth=2.0, label="Scaler $p_0$ fit", zorder=1)
    ax_bot.axhline(scaler_elec_current_p0*1.01, linestyle="--", color="darkred", linewidth=1.5, label="Scaler $p_0$ + 1%")
    ax_bot.axhline(scaler_elec_current_p0*0.99, linestyle="--", color="darkred", linewidth=1.5, label="Scaler $p_0$ - 1%")
    
    ax_bot.errorbar(scaler_elec_df["ibeam"].to_numpy(), scaler_elec_df["scaler_yield"].to_numpy(), yerr=scaler_elec_df["scaler_err"].to_numpy(), fmt="o", color="red", label="Scaler", zorder=2)
    ax_bot.set_ylabel("Electron Scaler Yield")
    ax_bot.set_xlabel("Beam Current (µA)")
    ax_bot.text(0.02, 0.95, f"$\\chi^2/ndf$ = {scaler_elec_chi2_ndf: .2f}", transform=ax_bot.transAxes, ha="left", va="top")
    ax_bot.legend(loc="upper left", bbox_to_anchor=(1.02, 1))

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    # -------------------------------------------------
    if have_positrons:
        fig, (ax_top, ax_bot) = plt.subplots(2, 1,gridspec_kw={"height_ratios": [1, 1], "hspace": 0.01},sharex=True)

        fig.suptitle(f"{selected_run_type.upper()} {selected_beam_pass}Pass {target_longname} Positron Yield vs Run Number", y=0.98)

        ax_top.axhline(pos_p0, linestyle="-", color="cornflowerblue", linewidth=2.0, label="Tracked $p_0$ fit", zorder=1)
        ax_top.axhline(pos_p0*1.01, linestyle="--", color="cornflowerblue", linewidth=1.5, label="Tracked $p_0$ + 1%", zorder=1)
        ax_top.axhline(pos_p0*0.99, linestyle="--", color="cornflowerblue", linewidth=1.5, label="Tracked $p_0$ - 1%", zorder=1)
        ax_top.errorbar(x_idx_pos, pos_df["yield"].to_numpy(), yerr=pos_df["yield_err"].to_numpy(), fmt="o",  color="red", label="Tracked", zorder=2)
        ax_top.set_ylabel("Positron Tracked Yield")
        ax_top.text(0.02, 0.95, f"$\\chi^2/ndf$ = {pos_chi2_ndf: .2f}", transform=ax_top.transAxes, ha="left", va="top")
        ax_top.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
        ax_top.tick_params(labelbottom=False)

        ax_bot.axhline(scaler_pos_p0, linestyle="-", color="darkred", linewidth=2.0, label="Scaler $p_0$ fit", zorder=1)
        ax_bot.axhline(scaler_pos_p0*1.01, linestyle="--", color="darkred", linewidth=1.5, label="Scaler $p_0$ + 1%", zorder=1)
        ax_bot.axhline(scaler_pos_p0*0.99, linestyle="--", color="darkred", linewidth=1.5, label="Scaler $p_0$ fit - 1%", zorder=1)
        ax_bot.errorbar(scaler_x_idx_pos, scaler_pos_df["scaler_yield"].to_numpy(), yerr=scaler_pos_df["scaler_err"].to_numpy(), fmt="o", color="orange", label="Scaler", zorder=2)
        ax_bot.set_xticks(x_idx_pos)
        ax_bot.set_xticklabels(pos_df["runnum"].astype(str), rotation=45)
        ax_bot.set_ylabel("Positron Scaler Yield")
        ax_bot.set_xlabel("Run Number")
        ax_bot.text(0.02, 0.95, f"$\\chi^2/ndf$ = {scaler_pos_chi2_ndf: .2f}", transform=ax_bot.transAxes, ha="left", va="top")
        ax_bot.legend(loc="upper left", bbox_to_anchor=(1.02, 1))

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # -------------------------------------------------
        fig, (ax_top, ax_bot) = plt.subplots(2, 1,gridspec_kw={"height_ratios": [1, 1], "hspace": 0.01},sharex=True)

        fig.suptitle(f"{selected_run_type.upper()} {selected_beam_pass}Pass {target_longname} Positron Yield vs Beam Current", y=0.98)

        ax_top.axhline(pos_current_p0, linestyle="-", color="cornflowerblue", linewidth=2.0, label="Tracked $p_0$ fit", zorder=1)
        ax_top.axhline(pos_current_p0*1.01, linestyle="--", color="cornflowerblue", linewidth=1.5, label="Tracked $p_0$ + 1%", zorder=1)
        ax_top.axhline(pos_current_p0*0.99, linestyle="--", color="cornflowerblue", linewidth=1.5, label="Tracked $p_0$ - 1%", zorder=1)
        
        ax_top.errorbar(pos_df["ibeam"].to_numpy(), pos_df["yield"].to_numpy(), yerr=pos_df["yield_err"].to_numpy(), fmt="o",  color="red", label="Tracked", zorder=2)
        ax_top.set_ylabel("Positron Tracked Yield")
        ax_top.text(0.02, 0.95, f"$\\chi^2/ndf$ = {pos_current_chi2_ndf: .2f}", transform=ax_top.transAxes, ha="left", va="top")
        ax_top.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
        ax_top.tick_params(labelbottom=False)

        ax_bot.axhline(scaler_pos_current_p0, linestyle="-", color="darkred", linewidth=2.0, label="Scaler $p_0$ fit", zorder=1)
        ax_bot.axhline(scaler_pos_current_p0*1.01, linestyle="--", color="darkred", linewidth=1.5, label="Scaler $p_0$ + 1%", zorder=1)
        ax_bot.axhline(scaler_pos_current_p0*0.99, linestyle="--", color="darkred", linewidth=1.5, label="Scaler $p_0$ fit", zorder=1)
        
        ax_bot.errorbar(scaler_pos_df["ibeam"].to_numpy(), scaler_pos_df["scaler_yield"].to_numpy(), yerr=scaler_pos_df["scaler_err"].to_numpy(), fmt="o", color="orange", label="Scaler", zorder=2)
        ax_bot.set_ylabel("Positron Scaler Yield")
        ax_bot.set_xlabel("Beam Current (µA)")
        ax_bot.text(0.02, 0.95, f"$\\chi^2/ndf$ = {scaler_pos_chi2_ndf: .2f}", transform=ax_bot.transAxes, ha="left", va="top")
        ax_bot.legend(loc="upper left", bbox_to_anchor=(1.02, 1))

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

print(f"Saved PDF → {pdf_filepath}")
