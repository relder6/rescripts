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
from scipy.optimize import curve_fit
from collections import defaultdict
from INIT.config import parse_run_type, parse_beam_pass, parse_target

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
USING_SYST_ERR_EST = True

arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None

selected_run_type = parse_run_type(arg1)
target_abbrev, target_longname, target_shortname, target_A, target_Z = parse_target(arg2)

os.makedirs("CSVs", exist_ok=True)
bc_csv = f"CSVs/{selected_run_type.upper()}_bin_centered_{target_abbrev}.csv"

os.makedirs("PDFs", exist_ok=True)
pdf_output = f"PDFs/{selected_run_type.upper()}_rosenbluth_separation_{target_abbrev}.pdf"
pp = PdfPages(pdf_output)

# -----------------------------------------------------
# Reading the csv
# -----------------------------------------------------
df = pd.read_csv(bc_csv)

required_cols = ["bin_num", "epsilon", "sigma_R", "sigma_R_err",
                 "bc_epsilon", "bc_sigma_R", "bc_sigma_R_err",
                 "xbj", "bc_xbj", "q2", "bc_q2"]

missing = [c for c in required_cols if c not in df.columns]

if missing:
    raise ValueError(f"Missing columns in {bc_csv}: {missing}")

# -----------------------------------------------------
# Readinging in bc_csv in dictionary by bin num
# -----------------------------------------------------
bc_data = {}
data = {}

for bin_num, sub in df.groupby("bin_num"):
    
    bc_data[bin_num] = {"epsilon": sub["bc_epsilon"].to_numpy(),
                        "sigma_R": sub["bc_sigma_R"].to_numpy(),
                        "sigma_R_err": sub["bc_sigma_R_err"].to_numpy(),
                        "xbj": sub["bc_xbj"].to_numpy(),
                        "q2": sub["bc_q2"].to_numpy()}
    
    data[bin_num] = {"epsilon": sub["epsilon"].to_numpy(),
                     "sigma_R": sub["sigma_R"].to_numpy(),
                     "sigma_R_err": sub["sigma_R_err"].to_numpy(),
                     "xbj": sub["xbj"].to_numpy(),
                     "q2": sub["q2"].to_numpy(),}

def linear_fit(x, intercept, slope):
    return intercept + slope * x

fit_results = []
eps_rows = []

with PdfPages(pdf_output) as pp:
    for bin_num in sorted(bc_data.keys()):

        eps_bc = np.array(bc_data[bin_num]["epsilon"])
        sig_bc = np.array(bc_data[bin_num]["sigma_R"])
        sig_err_bc = np.array(bc_data[bin_num]["sigma_R_err"])
        eps_raw = np.array(data[bin_num]["epsilon"])
        sig_raw = np.array(data[bin_num]["sigma_R"])
        sig_err_raw = np.array(data[bin_num]["sigma_R_err"])

        if USING_SYST_ERR_EST:
            syst_err = 0.015
        else:
            syst_err = 0
            
        sig_err_bc = np.sqrt(sig_err_bc**2 + (syst_err * sig_bc)**2)
        sig_err_raw = np.sqrt(sig_err_raw**2 + (syst_err * sig_raw)**2)
        
        xbj_arr = np.array(bc_data[bin_num]["xbj"])
        q2_arr = np.array(bc_data[bin_num]["q2"])
        xbjavg = np.mean(xbj_arr)
        n_points = len(q2_arr)
        q2avg = np.mean(q2_arr)
        q2max = np.max(q2_arr)
        q2min = np.min(q2_arr)
        q2_err = np.std(q2_arr, ddof=1) if n_points > 2 else (q2max - q2min) / 2

        eps_unique = defaultdict(lambda: {"sigR": [], "sigR_err": []})
        
        for eps, sig, sig_err in zip(eps_bc, sig_bc, sig_err_bc):
            eps_unique[eps]["sigR"].append(sig)
            eps_unique[eps]["sigR_err"].append(sig_err)

        eps_comb, sig_comb, sig_err_comb = [], [], []
        
        for eps, vals in eps_unique.items():
            sig_arr = np.array(vals["sigR"])
            sig_err_arr = np.array(vals["sigR_err"])
            w = 1.0 / sig_err_arr**2
            sig_avg = np.sum(sig_arr * w) / np.sum(w)
            sig_err = np.sqrt(1.0 / np.sum(w))
            eps_comb.append(eps)
            sig_comb.append(sig_avg)
            sig_err_comb.append(sig_err)

        eps_comb = np.array(eps_comb)
        sig_comb = np.array(sig_comb)
        sig_err_comb = np.array(sig_err_comb)

        popt, pcov = curve_fit(linear_fit, eps_comb, sig_comb, sigma=sig_err_comb, absolute_sigma=True)
        intercept, slope = popt
        int_err = np.sqrt(pcov[0, 0])
        slope_err = np.sqrt(pcov[1, 1])
        R = slope / intercept
        R_err = np.sqrt((pcov[1, 1] / intercept**2) + (slope**2 * pcov[0, 0] / intercept**4) - (2 * slope * pcov[0, 1] / intercept**3))

        fit_results.append({"bin_num": bin_num,
                            "xbj": xbjavg,
                            "q2_avg": q2avg,
                            "q2_err": q2_err,
                            "sigma_T": intercept,
                            "sigma_T_err": int_err,
                            "sigma_L": slope,
                            "sigma_L_err": slope_err,
                            "R": R,
                            "R_err": R_err,
                            "n_eps": len(eps_comb)})

        fig, ax = plt.subplots(figsize=(6, 5))

        ax.errorbar(eps_raw, sig_raw, yerr=sig_err_raw, fmt="o", color="gray", alpha = 0.4, markersize=3, label="pre bin-centering")

        ax.errorbar(eps_bc, sig_bc, yerr=sig_err_bc,fmt="o", color="navy", alpha = 0.8, markersize=3,label="bin-centered points")

        ax.errorbar(eps_comb, sig_comb, yerr = sig_err_comb, fmt = "o", color = "red", markersize = 5,label = "weighted avg per ε")

        eps_min = min(eps_raw.min(), eps_bc.min(), eps_comb.min())
        eps_max = max(eps_raw.max(), eps_bc.max(), eps_comb.max())
        x_fit = np.linspace(eps_min, eps_max, 200)
        y_fit = linear_fit(x_fit, intercept, slope)
        ax.plot(x_fit, y_fit, "--", color="red",label=("Linear fit\n"
                                                       rf"$\sigma_L = {slope:.4f} \pm {slope_err:.4f}$" "\n"
                                                       rf"$\sigma_T = {intercept:.4f} \pm {int_err:.4f}$" "\n"
                                                       rf"$R = {R:.4f} \pm {R_err:.4f}$""\n"f"(includes est. {syst_err} syst. err.)"))

        all_sig = np.concatenate([sig_raw, sig_bc, sig_comb])
        y_min, y_max = all_sig.min(), all_sig.max()
        margin = 0.1
        dy = margin * (y_max - y_min if y_max > y_min else 1.0)
        ax.set_xlim(eps_min - 0.05, eps_max + 0.05)
        ax.set_ylim(y_min - dy, y_max + dy)

        ax.set_title(f"{target_longname} Rosenbluth Separation\n"r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f} ± {q2_err:.3f}",fontsize=10)
        ax.set_xlabel(r"$\epsilon$")
        ax.set_ylabel(r"$d \sigma / d \Omega / dE' / \Gamma$ ($\mu$b/sr)")
        ax.legend(fontsize=8)
        ax.grid(True)

        pp.savefig(fig)
        plt.close(fig)

print(f"PDF saved to {pdf_output}")

csv_output = f"CSVs/{selected_run_type.upper()}_rosenbluth_fit_{target_abbrev}.csv"
pd.DataFrame(fit_results).to_csv(csv_output, index=False)
print(f"CSV of fits saved to {csv_output}")
