#!/usr/bin/env python3

import pandas as pd
import os, re, sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

input_directory = "../XSEC/DATA_to_MC"

output_csv = "CSVs/RME_ytars"

all_targets = ["Al", "C", "Cu", "LD2", "LH2"]

# targets = ["Al", "C", "Cu"]

beam_passes = ["4pass", "5pass"]

all_rows = []

for target in all_targets:
    for beam_pass in beam_passes:
        csv_file = f"{input_directory}/{target.upper()}/DATA_to_MC_hmsdis_{beam_pass}_{target.lower()}_H_gtr_y.csv"
        if not os.path.exists(csv_file):
            print(f"Skipping missing file: {csv_file}")
            continue

        df = pd.read_csv(csv_file)
        df["target"] = target
        df["beam_pass"] = beam_pass
        df["ytar"] = df["bin_center"]
        df["yield"] = df["yield_sub"]
        df["yield_err"] = df["err_sub"]
        df["mc"] = df["yield_mc"]
        df["mc_err"] = df["err_mc"]

        df_out = df[["target", "beam_pass", "ytar", "yield", "yield_err", "mc", "mc_err"]]

        df_out = df_out.dropna()

        all_rows.append(df_out)

df_all = pd.concat(all_rows, ignore_index = True)

df_all.to_csv(output_csv, index=False)

print(f"Saved compiled CSV → {output_csv}")

# ---------------------------------------------
# Plotting stuff
# ---------------------------------------------
def fit_peak(x, y, yerr, window=0.15):

    x = np.asarray(x)
    y = np.asarray(y)
    yerr = np.asarray(yerr)

    peak_idx = np.argmax(y)
    peak_x = x[peak_idx]
    peak_y = y[peak_idx]

    mask = (x >= peak_x - window) & (x <= peak_x + window)

    x_fit = x[mask]
    y_fit = y[mask]
    yerr_fit = yerr[mask]

    if x_fit.size < 3:
        return np.nan, np.nan, np.nan, True

    amp_guess = np.max(y_fit)
    mean_guess = peak_x
    sigma_guess = np.std(x_fit)

    def gaussian(x, amp, mean, sigma):
        return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

    try:
        popt, pcov = curve_fit(gaussian,x_fit, y_fit,p0=[amp_guess, mean_guess, sigma_guess],sigma=yerr_fit,absolute_sigma=True)
        amp_fit, mean_fit, sigma_fit = popt
    except Exception:
        return np.nan, np.nan, np.nan, True

    if not (-0.6 <= mean_fit <= 0.6):
        return np.nan, np.nan, np.nan, True

    return amp_fit, mean_fit, sigma_fit, False

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

for target in all_targets:
    for beam_pass in beam_passes:
        
        df_out = df_all[(df_all["target"] == target) & (df_all["beam_pass"] == beam_pass)].dropna()

        amp_y, mean_y, sigma_y, fail_y = fit_peak(df_out["ytar"].to_numpy(),
                                                  df_out["yield"].to_numpy(),
                                                  df_out["yield_err"].to_numpy(),)

        amp_mc, mean_mc, sigma_mc, fail_mc = fit_peak(df_out["ytar"].to_numpy(),
                                                      df_out["mc"].to_numpy(),
                                                      df_out["mc_err"].to_numpy(),)

        delta_peak = mean_y - mean_mc

        df_out["yield_amp_fit"] = amp_y
        df_out["yield_mean_fit"] = mean_y
        df_out["yield_sigma_fit"] = sigma_y

        df_out["mc_amp_fit"] = amp_mc
        df_out["mc_mean_fit"] = mean_mc
        df_out["mc_sigma_fit"] = sigma_mc

        df_out["yield_fit_failed"] = fail_y
        df_out["mc_fit_failed"] = fail_mc

        fig, ax = plt.subplots(figsize=(7, 6))

        x = df_out["ytar"].to_numpy()
        y = df_out["yield"].to_numpy()
        yerr = df_out["yield_err"].to_numpy()

        mc = df_out["mc"].to_numpy()
        mcerr = df_out["mc_err"].to_numpy()

        ax.errorbar(x, y, yerr=yerr, fmt="o", markersize=3, color="navy", label="Data")
        ax.errorbar(x, mc, yerr=mcerr, fmt="o", markersize=3,color="red", label="MC")

        window = 0.15

        x_fit_y = np.linspace(mean_y - window, mean_y + window, 200)
        x_fit_mc = np.linspace(mean_mc - window, mean_mc + window, 200)

        if not fail_y:
            ax.plot(x_fit_y,gaussian(x_fit_y, amp_y, mean_y, sigma_y),color="navy",lw=2)

        if not fail_mc:
            ax.plot(x_fit_mc,gaussian(x_fit_mc, amp_mc, mean_mc, sigma_mc),color="red",lw=2)

        ax.set_ylabel("Charge-Normalized Yields")
        ax.set_title(f"HMSDIS {beam_pass} {target} ytar Data vs MC Offset")
        ax.grid()
        ax.legend(loc="best",frameon=True, fontsize=9,title=f"$\Delta y = y_{{data}} - y_{{mc}}$ = {delta_peak:.4f}")

        plt.tight_layout()

        outfile = f"PNGs/ytar_offset_hmsdis_{beam_pass}_{target.lower()}.png"
        plt.savefig(outfile, dpi=200)

        plt.close()

        print(f"Saved plot → {outfile}")
