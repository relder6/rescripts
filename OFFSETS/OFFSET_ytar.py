#!/usr/bin/env python3

import pandas as pd
import os, re, sys
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

input_directory = "../XSEC/DATA_to_MC"
output_csv = "CSVs/ytar_offsets.csv"
output_png_dir = "PNGs"

solid_targets = ["Al", "C", "Cu", "Dummy_up", "Dummy_down"]
cryo_targets = ["LD2", "LH2"]
all_targets = solid_targets + cryo_targets

beam_passes = ["3pass", "4pass", "5pass"]
phases = ["phaseI", "phaseII"]

os.makedirs(os.path.dirname(output_csv), exist_ok=True)
os.makedirs(output_png_dir, exist_ok=True)

all_rows = []

def load_target_csv(target, beam_pass, phase):
    csv_file = f"{input_directory}/{target.upper()}/DATA_to_MC_hmsdis_{beam_pass}_{phase}_{target.lower()}_H_gtr_y.csv"
    if not os.path.exists(csv_file):
        print(f"Skipping missing file: {csv_file}")
        return None
    df = pd.read_csv(csv_file)
    df["target"] = target
    df["beam_pass"] = beam_pass
    df["phase"] = phase
    df["ytar"] = df["bin_center"]
    df["yield"] = df["yield_sub"]
    df["yield_err"] = df["err_sub"]
    df["mc"] = df["yield_mc"]
    df["mc_err"] = df["err_mc"]
    return df[["target", "beam_pass", "phase", "ytar", "yield", "yield_err", "mc", "mc_err"]].dropna()

def weighted_centroid(x, y, yerr, window=0.15):
    x = np.asarray(x)
    y = np.asarray(y)
    yerr = np.asarray(yerr)

    if x.size == 0:
        return np.nan

    peak_x = x[np.argmax(y)]
    m = (x >= peak_x - window) & (x <= peak_x + window) & np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr)
    if np.count_nonzero(m) < 3:
        return np.nan

    w = np.clip(y[m], 0, None) / np.maximum(yerr[m] ** 2, 1e-12)
    if not np.isfinite(w).any() or np.sum(w) <= 0:
        return np.nan

    return np.sum(x[m] * w) / np.sum(w)

def find_x_shift(x, y_data, y_mc, err_data, err_mc, bounds=(-1.0, 1.0)):
    x = np.asarray(x)
    y_data = np.asarray(y_data)
    y_mc = np.asarray(y_mc)
    err_data = np.asarray(err_data)
    err_mc = np.asarray(err_mc)

    mc_interp = interp1d(x, y_mc, kind="linear", bounds_error=False, fill_value="extrapolate")

    def chi2(dx):
        mc_shifted = mc_interp(x + dx)
        sigma2 = np.maximum(err_data**2 + err_mc**2, 1e-12)
        return np.sum((y_data - mc_shifted) ** 2 / sigma2)

    res = minimize_scalar(chi2, bounds=bounds, method="bounded")
    return res.x, res.fun

df_list = []
for target in all_targets:
    for beam_pass in beam_passes:
        for phase in phases:
            df_tmp = load_target_csv(target, beam_pass, phase)
            if df_tmp is None or df_tmp.empty:
                continue
            df_list.append(df_tmp)
if not df_list:
    raise RuntimeError("No input files were loaded.")

df_all = pd.concat(df_list, ignore_index = True)
df_all.to_csv(output_csv, index=False)
print(f"Saved compiled CSV → {output_csv}")

# ---------------------------------------------
# Plotting stuff
# ---------------------------------------------
for target in all_targets:
    for beam_pass in beam_passes:
        for phase in phases:
            df_out = df_all[(df_all["target"] == target) & (df_all["beam_pass"] == beam_pass) & (df_all["phase"] == phase)].dropna()
            if df_out.empty:
                continue

            x = df_out["ytar"].to_numpy()
            y = df_out["yield"].to_numpy()
            yerr = df_out["yield_err"].to_numpy()
            mc = df_out["mc"].to_numpy()
            mcerr = df_out["mc_err"].to_numpy()

            if target in cryo_targets:
                shift_value, chi2_min = find_x_shift(x, y, mc, yerr, mcerr, bounds=(-1.0, 1.0))
                mc_interp = interp1d(x, mc, kind="linear", bounds_error=False, fill_value="extrapolate")
                mc_shifted = mc_interp(x + shift_value)
                mc_plot = mc_shifted
                mc_plot_label = "MC (shifted)"
                mc_plot_color = "red"
            else:
                shift_data = weighted_centroid(x, y, yerr, window=0.10)
                shift_mc = weighted_centroid(x, mc, mcerr, window=0.10)
                if np.isfinite(shift_data) and np.isfinite(shift_mc):
                    shift_value = shift_data - shift_mc
                else:
                    shift_value, chi2_min = find_x_shift(x, y, mc, yerr, mcerr, bounds=(-0.30, 0.30))
                shift_label = f"Δy = {shift_value:.4f}"
                mc_plot = mc
                mc_plot_label = "MC"
                mc_plot_color = "red"

            fig, ax = plt.subplots(figsize=(7, 6))

            shift_label = f"Data-MC = {shift_value:.4f}"

            mc_unshifted = mc

            ax.errorbar(x, y, yerr=yerr, fmt="o", markersize=3, color="navy", label="Data")
            ax.errorbar(x, mc_plot, yerr=mcerr, fmt="o", markersize=3, color=mc_plot_color, label=mc_plot_label)
            ax.errorbar(x, mc_unshifted, yerr=mcerr, fmt="o", markersize=3, color="lightgrey", label="MC (unshifted)")

            if target not in cryo_targets:
                data_centroid = weighted_centroid(x, y, yerr, window=0.10)
                mc_centroid = weighted_centroid(x, mc_unshifted, mcerr, window=0.10)

                if np.isfinite(data_centroid):
                    ax.axvline(data_centroid, linestyle="--", color="navy", alpha=0.8, label=f"Data centroid = {data_centroid:.4f}")
                if np.isfinite(mc_centroid):
                    ax.axvline(mc_centroid, linestyle="--", color="lightgrey", alpha=0.9, label=f"MC centroid = {mc_centroid:.4f}")

            ax.set_ylabel("Charge-Normalized Yields")
            ax.set_title(f"{beam_pass} {target} ytar Data vs MC")
            ax.grid(True)
            ax.legend(loc="best", frameon=True, fontsize=9, title=shift_label)
            plt.tight_layout()

            outfile = f"{output_png_dir}/ytar_offset_hmsdis_{beam_pass}_{phase}_{target.lower()}.png"
            plt.savefig(outfile, dpi=200)
            plt.close()
            print(f"Saved plot -> {outfile}")
