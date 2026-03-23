#!/usr/bin/env python3

import pandas as pd
import os, re, sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# -----------------------------------------------------
# Input Settings and Such
# -----------------------------------------------------
ratio_directory = "../XSEC/DATA_to_MC"

beam_passes = {"4pass", "5pass"}

# targets = {"C", "Cu", "LD2", "LH2"}

targets = {"C", "Cu"}

# targets = {"LD2", "LH2"}

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

        df["ratio_offset_err"] = np.sqrt(df["ratio_err"]**2 + df["p0_fit_err"]**2)

        # df["ratio_offset_err"] = df["ratio_err"]

        df_out = df[["target", "beam_pass", "delta", "ratio", "ratio_err", "p0_fit", "ratio_offset", "ratio_offset_err"]]

        df_out = df_out.dropna()

        all_rows.append(df_out)

df_all = pd.concat(all_rows, ignore_index = True)

output_csv = "DELTA_corr_testing.csv"

df_all.to_csv(output_csv, index = False)

print(f"Saved compiled CSV → {output_csv}")

# -----------------------------------------------------
# 
# -----------------------------------------------------
x = df_all["delta"].values
y = df_all["ratio_offset"].values
sigma = df_all["ratio_offset_err"].values

mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
x_fit = x[mask]
y_fit = y[mask]
sigma_fit = sigma[mask]

dx = 1.0 # step in delta

s0 = np.sum(dx * np.ones_like(x_fit))
s1 = np.sum(x_fit * dx)
s2 = np.sum(x_fit**2 * dx)
s3 = np.sum(x_fit**3 * dx)
s4 = np.sum(x_fit**4 * dx)

total_integral = np.sum(y_fit * dx)

def poly4_fit_constr(x, b, c, d, e):
    a = (total_integral - b*s1 - c*s2 -d*s3 -e*s4) / s0
    return a + b * x + c * x**2 + d * x**3 + e * x**4

try:
    popt_free, pcov = curve_fit(poly4_fit_constr, x_fit, y_fit, sigma=sigma_fit, absolute_sigma = True)
    b, c, d, e = popt_free
    a = (total_integral - b*s1 - c*s2 - d*s3 - e*s4) / s0
    popt = np.array([a, b, c, d, e])
    p4_err = np.sqrt(np.diag(pcov))
except Exception as e:
    print(f"Fit failed for {csv_file}: {e}")
    popt_free = np.full(4, np.nan)
    popt = np.full(5, np.nan)
    p4_err = np.full(5, np.nan)

y_model = poly4_fit_constr(x_fit, *popt_free)
chi2 = np.sum(((y_fit - y_model)/sigma_fit)**2)
ndf = len(x_fit) - len(popt_free)
chi2_ndf = chi2 / ndf

x_smooth = np.linspace(np.min(x_fit), np.max(x_fit), 200)
y_smooth = poly4_fit_constr(x_smooth, *popt_free)

coeff_labels = ["a","b","c","d","e"]
coeff_text = f"p4 fit: $a + bx + cx^2 + dx^3 + ex^4$\n"
for label, val in zip(coeff_labels, popt):
    coeff_text += f"{label} = {val:.6e}\n"
coeff_text += f"Chi2 / ndf = {chi2_ndf:.2f} ({chi2:.2f}/{ndf})"

print(f"{coeff_text}")

integral_model = np.sum(poly4_fit_constr(x_fit, *popt_free) * dx)
print("Idata =", total_integral, "Imodel =", integral_model, "Imodel/Idata =", integral_model / total_integral)

# -----------------------------------------------------
# Plot data and fit
# -----------------------------------------------------
target_colors = {"C": "blue", "LD2": "green", "Cu": "red", "LH2": "orange"}

beam_markers = {"4pass": "o", "5pass": "x"}

for target in df_all["target"].unique():
    for beam in df_all["beam_pass"].unique():
        df_plot = df_all[(df_all["target"]==target) & (df_all["beam_pass"]==beam)]
        x = df_plot["delta"].to_numpy()
        y = df_plot["ratio_offset"].to_numpy()
        yerr = df_plot["ratio_offset_err"].to_numpy()
        plt.errorbar(x, y, yerr=yerr,fmt=beam_markers[beam],color=target_colors[target],label=f"{target} {beam}",alpha=0.7)
plt.text(0.05, 0.95, coeff_text,transform=plt.gca().transAxes,fontsize=10,verticalalignment='top',bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.8))
plt.plot(x_smooth, y_smooth, "r--", label="p(4) fit", linewidth=2)
plt.xlabel("Delta", fontsize = 14)
plt.ylabel("Data/MC Ratio", fontsize = 14)
plt.title("Delta Correction Studies", fontsize = 16)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


