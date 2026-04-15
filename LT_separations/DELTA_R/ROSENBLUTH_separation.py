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
from INIT.config import get_data_cuts, get_common_values
from scipy.optimize import curve_fit
from collections import defaultdict
from INIT.config import parse_run_type, parse_beam_pass, parse_target, parse_bins

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
USING_SYST_ERR_EST = True

arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None
arg3 = sys.argv[3] if len(sys.argv) > 3 else None

selected_run_type = parse_run_type(arg1)
num_abbrev, num_longname, num_shortname, num_A, num_Z = parse_target(arg2)
denom_abbrev, denom_longname, denom_shortname, denom_A, denom_Z = parse_target(arg3)

A_ratio = num_A / denom_A
num_N = num_A - num_Z
denom_N = denom_A - denom_Z

vals = get_common_values()
ebeam_4pass = vals["ebeam_4pass"]
theta_4pass = vals["angle_4pass"]
ebeam_5pass = vals["ebeam_5pass"]
theta_5pass = vals["angle_5pass"]

model_xsec_dir = "../../../mc-single-arm/util/dis_xec"

xsec_ratio_dir = "../../XSEC/FORM_xsec/RATIOS"

beam_passes = ["4pass", "5pass"]

os.makedirs("CSVs", exist_ok=True)
bc_csv = f"CSVs/DELTA_R_{selected_run_type.upper()}_bin_centered_{num_abbrev}_to_{denom_abbrev}.csv"

os.makedirs("PDFs", exist_ok=True)
pdf_output = f"PDFs/DELTA_R_{selected_run_type.upper()}_rosenbluth_separation_{num_abbrev}_to_{denom_abbrev}.pdf"
pp = PdfPages(pdf_output)

# -----------------------------------------------------
# Reading the csv
# -----------------------------------------------------
df = pd.read_csv(bc_csv)

required_cols = ["bin_num", "xbj", "q2", "epsilon_p", "sigma_num_to_sigma_denom", "sigma_num_to_sigma_denom_err",
                 "bc_xbj", "bc_q2", "bc_epsilon_p", "bc_sigma_num_to_sigma_denom", "bc_sigma_num_to_sigma_denom_err",]

missing = [c for c in required_cols if c not in df.columns]

if missing:
    raise ValueError(f"Missing columns in {bc_csv}: {missing}")

# -----------------------------------------------------
# Readinging in bc_csv in dictionary by bin num
# -----------------------------------------------------
bc_data = {}
data = {}

for bin_num, sub in df.groupby("bin_num"):
    
    bc_data[bin_num] = {"epsilon_p": sub["bc_epsilon_p"].to_numpy(),
                        "xsec_ratio": sub["bc_sigma_num_to_sigma_denom"].to_numpy(),
                        "xsec_ratio_err": sub["bc_sigma_num_to_sigma_denom_err"].to_numpy(),
                        "xbj": sub["bc_xbj"].to_numpy(),
                        "q2": sub["bc_q2"].to_numpy(),}
    
    data[bin_num] = {"epsilon_p": sub["epsilon_p"].to_numpy(),
                     "xsec_ratio": sub["sigma_num_to_sigma_denom"].to_numpy(),
                     "xsec_ratio_err": sub["sigma_num_to_sigma_denom_err"].to_numpy(),
                     "xbj": sub["xbj"].to_numpy(),
                     "q2": sub["q2"].to_numpy(),}

# def linear_fit(x, intercept, slope):
#     return intercept + slope * x

def modified_linear_fit(x, sig_t_ratio, deltaR):
    return sig_t_ratio * (1 + ( deltaR * x))

fit_results = []
epsp_rows = []

with PdfPages(pdf_output) as pp:
    for bin_num in sorted(bc_data.keys()):

        epsp_bc = np.array(bc_data[bin_num]["epsilon_p"])
        ratio_bc = np.array(bc_data[bin_num]["xsec_ratio"])
        ratio_err_bc = np.array(bc_data[bin_num]["xsec_ratio_err"])
        epsp_raw = np.array(data[bin_num]["epsilon_p"])
        ratio_raw = np.array(data[bin_num]["xsec_ratio"])
        ratio_err_raw = np.array(data[bin_num]["xsec_ratio_err"])

        if USING_SYST_ERR_EST:
            syst_err = 0.010
        else:
            syst_err = 0
            
        ratio_err_bc = np.sqrt(ratio_err_bc**2 + (syst_err * ratio_bc)**2)
        ratio_err_raw = np.sqrt(ratio_err_raw**2 + (syst_err * ratio_raw)**2)

        xbj_arr = np.array(bc_data[bin_num]["xbj"])
        q2_arr = np.array(bc_data[bin_num]["q2"])
        xbjavg = np.mean(xbj_arr)
        n_points = len(q2_arr)
        q2avg = np.mean(q2_arr)
        q2max = np.max(q2_arr)
        q2min = np.min(q2_arr)
        q2_err = np.std(q2_arr, ddof=1) if n_points > 2 else (q2max - q2min) / 2

        epsp_unique = defaultdict(lambda: {"ratio": [], "ratio_err": []})

        for epsp, rat, rat_err in zip(epsp_bc, ratio_bc, ratio_err_bc):
            epsp_unique[epsp]["ratio"].append(rat)
            epsp_unique[epsp]["ratio_err"].append(rat_err)

        epsp_comb, rat_comb, rat_err_comb = [], [], []

        for epsp, vals in epsp_unique.items():
            rat_arr = np.array(vals["ratio"])
            rat_err_arr = np.array(vals["ratio_err"])
            w = 1.0 / rat_err_arr**2
            rat_avg = np.sum(rat_arr * w) / np.sum(w)
            rat_err = np.sqrt(1.0 / np.sum(w))
            epsp_comb.append(epsp)
            rat_comb.append(rat_avg)
            rat_err_comb.append(rat_err)

        epsp_comb = np.array(epsp_comb)
        rat_comb = np.array(rat_comb)
        rat_err_comb = np.array(rat_err_comb)

        popt, pcov = curve_fit(modified_linear_fit, epsp_comb, rat_comb, sigma=rat_err_comb, absolute_sigma=True)
        sig_t_ratio, deltaR = popt
        sig_t_ratio_err = np.sqrt(pcov[0, 0])
        deltaR_err = np.sqrt(pcov[1, 1])

        fit_results.append({"bin_num": bin_num,
                            "xbj": xbjavg,
                            "q2_avg": q2avg,
                            "q2_err": q2_err,
                            "delta_R": deltaR,
                            "delta_R_err": deltaR_err,
                            "sigma_t_ratio": sig_t_ratio,
                            "sigma_t_ratio_err": sig_t_ratio_err})

        fig, ax = plt.subplots(figsize = (6, 5))

        ax.errorbar(epsp_raw, ratio_raw, yerr = ratio_err_raw, fmt = "o", color = "gray", alpha = 0.4, markersize = 3, label = "pre bin-centering")

        ax.errorbar(epsp_bc, ratio_bc, yerr = ratio_err_bc, fmt = "o", color = "navy", alpha = 0.8, markersize = 3,label = "bin-centered points")

        ax.errorbar(epsp_comb, rat_comb, yerr = rat_err_comb, fmt = "o", color = "red", markersize = 5,label = "weighted avg per ε'")

        epsp_min = min(epsp_raw.min(), epsp_bc.min(), epsp_comb.min())
        epsp_max = max(epsp_raw.max(), epsp_bc.max(), epsp_comb.max())
        x_fit = np.linspace(epsp_min, epsp_max, 200)
        y_fit = modified_linear_fit(x_fit, sig_t_ratio, deltaR)
        ax.plot(x_fit, y_fit, "--", color="red",label=(r"Modified Linear Fit: $\frac{\sigma_A}{\sigma_D} = \frac{\sigma_A^T}{\sigma_D^T} \left(1 + \Delta R \epsilon' \right)$"f"\n"f"(includes est. {syst_err} syst. err.)"))
        
        all_rat = np.concatenate([ratio_raw, ratio_bc, rat_comb])
                
        y_min, y_max = all_rat.min(), all_rat.max()
                
        margin = 0.1
        dy = margin * (y_max - y_min if y_max > y_min else 1.0)
        ax.set_xlim(epsp_min - 0.05, epsp_max + 0.05)
        ax.set_ylim(y_min - dy, y_max + dy)

        ax.set_title(f"{num_longname}/{denom_longname} Rosenbluth Separation\n"r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f} ± {q2_err:.3f}\n $\Delta$R = {deltaR:.4f} ± {deltaR_err:.4f},"f"\t"r"$\sigma_A^T \backslash \sigma_D^T =$"f"{sig_t_ratio:.4f}"r" $\pm$"f"{sig_t_ratio_err:.4f}", fontsize=10)
        ax.set_xlabel(r"$\epsilon$'")
        ax.set_ylabel(r"$ \sigma_A / \sigma_D$")
        
        ax.legend(fontsize=8)
        ax.grid(True)

        pp.savefig(fig)
        plt.close(fig)

print(f"PDF saved to {pdf_output}")
csv_output = f"CSVs/DELTA_R_{selected_run_type.upper()}_rosenbluth_fit_{num_abbrev}_to_{denom_abbrev}.csv"
pd.DataFrame(fit_results).to_csv(csv_output, index=False)
print(f"CSV of fits saved to {csv_output}")
