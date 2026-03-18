#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import os, re, sys, csv
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import get_common_run_inputs, get_data_cuts, get_common_values
from scipy.optimize import curve_fit

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
if len(sys.argv) == 6:
    selected_run_type = sys.argv[1].strip().lower()
    selected_beam_pass = sys.argv[2].strip()
    selected_num = sys.argv[3].strip().lower()
    selected_denom = sys.argv[4].strip().lower()
    nbins = int(sys.argv[5])
else:
    selected_run_type = input("Enter desired run type (default HMSDIS): ").strip().lower()
    if not selected_run_type:
        selected_run_type = "hmsdis"
    selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ").strip()
    selected_num = input("Enter numerator target: ").strip().lower()
    selected_denom = input("Enter denominator target: ").strip().lower()
    nbins = int(input("Enter bin number: "))
    
selected_beam_pass_to_energy_prefix = {
    "1": "2.", "2": "4.", "3": "6.", "4": "8.", "5": "10."
}
beam_prefix = selected_beam_pass_to_energy_prefix.get(selected_beam_pass)
if not beam_prefix:
    print(f"Unknown pass: {selected_beam_pass}.  Please try again.")
    exit(1)

selected_target_shortcut_to_target_variable = {
    "al":"al","al13":"al","aluminum":"al",
    "c":"c","c12":"c","carbon":"c",
    "cu":"cu","cu29":"cu","copper":"cu",
    "opt1":"optics1","optics1":"optics1",
    "opt2":"optics2","optics2":"optics2",
    "d2":"ld2","ld2":"ld2",
    "h2":"lh2","lh2":"lh2",
    "hole":"hole","chole":"hole","c-hole":"hole",
    "dummy":"dummy","dum":"dummy"
}

selected_target_shortname_to_title_longname = {
    "al":"Aluminum",
    "c":"Carbon",
    "cu":"Copper",
    "opt1":"Optics1",
    "opt2":"Optics2",
    "ld2":"Deuterium",
    "lh2":"Hydrogen",
    "hole":"Carbon Hole",
    "dummy":"Dummy"}

# Targets
num_short = selected_target_shortcut_to_target_variable.get(selected_num)
if not num_short:
    print(f"Unknown target: {selected_num}. Please try again.")
    exit(1)
num_long = selected_target_shortname_to_title_longname[num_short]
    
denom_short = selected_target_shortcut_to_target_variable.get(selected_denom)
if not denom_short:
    print(f"Unknown target: {selected_denom}. Please try again.")
    exit(1)
denom_long = selected_target_shortname_to_title_longname[denom_short]

selected_target_shortname_to_AZ = {"al":   (27, 13),
                                   "c":    (12, 6),
                                   "cu":   (64, 29),
                                   "ld2":  (2, 1),
                                   "lh2":  (1, 1)}

A_num, Z_num = selected_target_shortname_to_AZ.get(num_short)

N_num = A_num - Z_num

A_denom, Z_denom = selected_target_shortname_to_AZ.get(denom_short, (0.0, 0.0))

A_ratio = A_num / A_denom

vals = get_common_values()
ebeam_4pass = vals["ebeam_4pass"]
theta_4pass = vals["angle_4pass"]
ebeam_5pass = vals["ebeam_5pass"]
theta_5pass = vals["angle_5pass"]

model_xsec_dir = "../../../mc-single-arm/util/dis_xec"

xsec_ratio_dir = "../../XSEC/FORM_xsec/RATIOS"

beam_passes = ["4pass", "5pass"]

# -----------------------------------------------------
# Files
# -----------------------------------------------------
os.makedirs("CSVs", exist_ok=True)
bc_csv = f"CSVs/DELTA_R_{selected_run_type.upper()}_bin_centered_{num_short}_to_{denom_short}.csv"

os.makedirs("PDFs", exist_ok=True)
pdf_output = f"PDFs/DELTA_R_{selected_run_type.upper()}_rosenbluth_separation_{num_short}_to_{denom_short}.pdf"
pp = PdfPages(pdf_output)

# -----------------------------------------------------
# Readinging in bc_csv in dictionary by bin num
# -----------------------------------------------------
bin_data = {}

with open(bc_csv, "r", newline="") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        try:
            bin_num = int(row["bin_num"])
            epsilon_p = float(row["bc_epsilon_p"])
            xsec_ratio = float(row["bc_sigma_num_to_sigma_denom"])
            xsec_ratio_err = float(row["bc_sigma_num_to_sigma_denom_err"])
            xbj = float(row["bc_xbj"])
            q2 = float(row["bc_q2"])

            if bin_num not in bin_data:
                bin_data[bin_num] = {"epsilon_p": [],
                                     "xsec_ratio": [],
                                     "xsec_ratio_err": [],
                                     "xbj": [],
                                     "q2": [],}
            bin_data[bin_num]["epsilon_p"].append(epsilon_p)
            bin_data[bin_num]["xsec_ratio"].append(xsec_ratio)
            bin_data[bin_num]["xsec_ratio_err"].append(xsec_ratio_err)
            bin_data[bin_num]["xbj"].append(xbj)
            bin_data[bin_num]["q2"].append(q2)
            
        except KeyError:
            continue
        except ValueError:
            continue

def linear_fit(x, intercept, slope):
    return intercept + slope * x

fit_results = []

with PdfPages(pdf_output) as pp:
    print(f"Found {len(bin_data)} bins. Processing...")

    for bin_num, data in bin_data.items():
        epsp = np.array(data["epsilon_p"])
        ratio = np.array(data["xsec_ratio"])
        ratio_err = np.array(data["xsec_ratio_err"])

        epsp_unique = {}
        for e, r, rerr in zip(epsp, ratio, ratio_err):
            if e not in epsp_unique:
                epsp_unique[e] = {"ratio": [], "ratio_err": []}
            epsp_unique[e]["ratio"].append(r)
            epsp_unique[e]["ratio_err"].append(rerr)

        epsp = []
        ratio = []
        ratio_err = []

        for e, vals in epsp_unique.items():
            r = np.array(vals["ratio"])
            rerr = np.array(vals["ratio_err"])

            weights = 1 / rerr**2
            r_avg = np.sum(r * weights) / np.sum(weights)
            r_err = np.sqrt(1 / np.sum(weights))

            epsp.append(e)
            ratio.append(r_avg)
            ratio_err.append(r_err)

        epsp = np.array(epsp)
        ratio = np.array(ratio)
        ratio_err = np.array(ratio_err)

        xbjavg = np.mean(data["xbj"])
        n_points = len(data["q2"])
        q2avg = np.mean(data["q2"])
        q2max = np.max(data["q2"])
        q2min = np.min(data["q2"])
        q2_err = np.std(data["q2"], ddof=1) if n_points > 2 else (q2max - q2min)/2

        popt, pcov = curve_fit(linear_fit, epsp, ratio, sigma=ratio_err, absolute_sigma=True)
        intercept, slope = popt
        int_err = np.sqrt(pcov[0,0])
        slope_err = np.sqrt(pcov[1,1])
        R = float(slope) / float(intercept)
        R_err = np.sqrt((pcov[1,1] / intercept**2) + (slope**2 * pcov[0,0] / intercept**4) - (2 * slope * pcov[0,1] / intercept**3))

        print(f"bin {bin_num}: xbj = {xbjavg:.3f}, q2 = {q2avg:.3f} ± {q2_err:.3f}")
        print(f"σ_A - σ_T = {slope:.4f} ± {slope_err:.4f}")

        fit_results.append({"bin_num": bin_num,
                            "xbj": xbjavg,
                            "q2_avg": q2avg,
                            "q2_err": q2_err,
                            "delta_R": slope,
                            "delta_R_err": slope_err})
        
        fig, ax = plt.subplots(figsize=(6, 5))
        
        ax.errorbar(epsp, ratio, yerr=ratio_err, fmt="o", color="navy", markersize = 4)

        # Fit evaluated across full x-axis [0, 1]
        x_full = np.linspace(0.0, 1.0, 200)
        y_full = linear_fit(x_full, intercept, slope)
        ax.plot(x_full,y_full,"--",color="red",label=("Linear fit: y = mx + b\n"
                                                      rf"$\Delta R = {slope:.4f} \pm {slope_err:.4f}$"))

        ax.set_title(f"{num_long}/{denom_long} Rosenbluth Separation\n"r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f} ± {q2_err:.3f}\n $\Delta$R = {slope:.4f} ± {slope_err:.4f}", fontsize=10)
        ax.set_xlabel(r"$\epsilon$'")
        ax.set_ylabel(r"$ \sigma_A / \sigma_D$")
        ax.set_xlim(0, 1)
        ax.legend(fontsize=8)
        ax.grid(True)

        pp.savefig(fig)
        plt.close(fig)

print(f"PDF saved to {pdf_output}")
csv_output = f"CSVs/DELTA_R_{selected_run_type.upper()}_rosenbluth_fit_{num_short}_to_{denom_short}.csv"
pd.DataFrame(fit_results).to_csv(csv_output, index=False)
print(f"CSV of fits saved to {csv_output}")
