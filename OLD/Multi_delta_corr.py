#!/usr/bin/env python3

import pandas as pd
import os, re, sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# -----------------------------------------------------
# Input Settings and Such
# -----------------------------------------------------
orders = [3, 4, 5, 6]

ratio_directory = "../XSEC/DATA_to_MC"

beam_passes = {"4pass", "5pass"}

# beam_passes = {"5pass"}

# beam_passes = {"5pass"}

# targets = {"C", "Cu", "LD2", "LH2"}

# targets = {"LD2", "LH2"}

targets = {"C"}

# -----------------------------------------------------
# Reading in and compiling the csv; p0 fit to offset and overlay
# -----------------------------------------------------
all_rows = []

def poly0_fit(x, p0):
    return p0

for beam_pass in beam_passes:
    for target in targets:
        csv_file = f"{ratio_directory}/{target.upper()}/DATA_to_MC_hmsdis_{beam_pass}_{target.lower()}_H_gtr_dp.csv"

        if not os.path.exists(csv_file):
            print(f"Skipping missing file: {csv_file}")
            continue

        df = pd.read_csv(csv_file)

        df["target"] = f"{target}"
        df["beam_pass"] = f"{beam_pass}"
        df["delta"] = df["bin_center"]
        
        x = df["delta"].values
        y = df["ratio"].values
        sigma = df["ratio_err"].values

        mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
        x_fit = x[mask]
        y_fit = y[mask]
        
        sigma_fit = sigma[mask]

        try:
            popt, pcov = curve_fit(poly0_fit, x_fit, y_fit, sigma=sigma_fit, absolute_sigma = True)
            p0_fit = popt[0]
            p0_err = np.sqrt(np.diag(pcov))[0]
        except Exception as e:
            print(f"Fit failed for {csv_file}: {e}")
            p0_fit = p0_err = np.nan

        df["p0_fit"] = p0_fit

        df["p0_fit_err"] = p0_err

        df["p0_offset"] = 1 - df["p0_fit"]

        df["ratio_offset"] = df["ratio"] + df["p0_offset"]

        # df["ratio_offset_err"] = np.sqrt(df["ratio_err"]**2 + df["p0_fit_err"]**2)

        df["ratio_offset_err"] = df["ratio_err"]

        df_out = df[["target", "beam_pass", "delta", "ratio", "ratio_err", "p0_fit", "ratio_offset", "ratio_offset_err"]]

        df_out = df_out.dropna()

        all_rows.append(df_out)

df_all = pd.concat(all_rows, ignore_index = True)

output_csv = "CSVs/delta_corr_testing.csv"

df_all.to_csv(output_csv, index = False)

print(f"Saved compiled CSV → {output_csv}")
print("*********************************")

# -----------------------------------------------------
# Doing the thing super cool
# -----------------------------------------------------
x = df_all["delta"].values
y = df_all["ratio_offset"].values
sigma = df_all["ratio_offset_err"].values

mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
x_fit = x[mask]
y_fit = y[mask]
sigma_fit = sigma[mask]

dx = 1.0 # step in delta

for order in orders:
    s = np.array([np.sum(x_fit**k * dx) for k in range(order+1)])
    total_integral = np.sum(y_fit * dx)

    def polyN_fit_constr(x, *coeffs):
        a = (total_integral - np.sum([c * s[k+1] for k, c in enumerate(coeffs)])) / s[0]
        return a + sum(c * x**(k+1) for k, c in enumerate(coeffs))

    if len(x_fit) == 0:
        raise RuntimeError("No valid data points after mask.")
    if order < 1:
        raise RuntimeError("Polynomial order must be at least 1.")

    nfree = order
    p0 = np.ones(order)
    popt_free = np.full(order, np.nan)  
    popt     = np.full(order+1, np.nan)
    pN_err   = np.full(order+1, np.nan)

    try:
        popt_free, pcov = curve_fit(polyN_fit_constr, x_fit, y_fit, p0=p0, sigma=sigma_fit, absolute_sigma = True)
        popt[0] = (total_integral - np.sum([p * s[k+1] for k, p in enumerate(popt_free)])) / s[0]
        popt[1:] = popt_free
        pN_err = np.sqrt(np.diag(pcov))

    except Exception as e:
        print(f"Fit failed for {csv_file}: {e}")

    y_model = polyN_fit_constr(x_fit, *popt_free)
    chi2 = np.sum(((y_fit - y_model)/sigma_fit)**2)
    ndf = len(x_fit) - len(popt_free)
    chi2_ndf = chi2 / ndf

    x_smooth = np.linspace(np.min(x_fit), np.max(x_fit), 200)
    y_smooth = polyN_fit_constr(x_smooth, *popt_free)

    coeff_names = ["a"] + [chr(ord("b")+i) for i in range(order)]
    
    coeff_text = (f"p{order} fit: $a + bx + cx^2 + \\dots + {coeff_names[-1]}x^{order}$\n")

    for label, val in zip(coeff_names, popt):
        coeff_text += f"{label} = {val:.6e}\n"

    coeff_text += f"Chi2 / ndf = {chi2_ndf:.2f} ({chi2:.2f}/{ndf})"

    print(f"{coeff_text}")

    integral_model = np.sum(polyN_fit_constr(x_fit, *popt_free) * dx)

    print("int data =", total_integral, "int model =", integral_model, "int_model/int_data =", integral_model / total_integral)

    print("")

    plt.figure(figsize=(10, 6))

    target_colors = {"C": "blue", "LD2": "green", "Cu": "red", "LH2": "orange"}
    
    beam_markers = {"4pass": "o", "5pass": "x"}

    for target in df_all["target"].unique():
        for beam in df_all["beam_pass"].unique():
            df_plot = df_all[(df_all["target"]==target) & (df_all["beam_pass"]==beam)]

            if df_plot.empty:
                continue
        
            x = df_plot["delta"].to_numpy()
            y = df_plot["ratio_offset"].to_numpy()
            yerr = df_plot["ratio_offset_err"].to_numpy()
        
            plt.errorbar(x, y, yerr=yerr,fmt=beam_markers[beam],color=target_colors[target],label=f"{target} {beam}",alpha=0.7)

    plt.plot(x_smooth, y_smooth, "r--", label=f"p{order} fit", linewidth=2)
        
    plt.text(0.05, 0.95, coeff_text,transform=plt.gca().transAxes,fontsize=10,verticalalignment='top',bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.8))

    plt.xlabel("Delta", fontsize=14)
    plt.ylabel("Data/MC Ratio", fontsize=14)
    plt.title(f"Delta Correction Studies: p{order} Fit", fontsize=16)
    # plt.title("DATA/MC Ratio Comparison After Delta Correction", fontsize=16)
    plt.xlim(-8.0, 8.0)
    plt.ylim(0.912, 1.068)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    residuals = (y_fit - y_model) / y_model
    residual_err = sigma_fit / y_model

    def p0_resid(x, c):
        return c

    try:
        popt_r, pcov_r = curve_fit(p0_resid,x_fit,residuals,sigma=residual_err,absolute_sigma=True)
        c_fit = popt_r[0]
        c_err = np.sqrt(np.diag(pcov_r))[0]
    except Exception as e:
        print("Residual fit failed:", e)
        c_fit, c_err = np.nan, np.nan

    plt.figure(figsize=(10,4))
    x_line = np.linspace(np.min(x_fit), np.max(x_fit), 100)
    plt.plot(x_line, np.full_like(x_line, c_fit), 'r--',
             label=f"p0 fit: {c_fit:.3f} ± {c_err:.3f}")
    plt.errorbar(x_fit, residuals, yerr=residual_err, fmt = "o", color = "navy", alpha=0.7)
    plt.xlabel("Delta")
    plt.ylabel("Residuals")
    plt.title(f"Residuals vs Delta (p{order})")
    plt.grid(True)
    plt.legend()

plt.show()
