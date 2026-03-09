#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import os, re, sys, csv
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import get_common_run_inputs
from scipy.optimize import curve_fit

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, selected_target_titlename, selected_target_A, selected_target_Z = get_common_run_inputs()

# -----------------------------------------------------
# Files
# -----------------------------------------------------
os.makedirs("CSVs", exist_ok=True)
bc_csv = f"CSVs/{selected_run_type.upper()}_bin_centered_{selected_target_shortname}.csv"

os.makedirs("PDFs", exist_ok=True)
pdf_output = f"PDFs/{selected_run_type.upper()}_rosenbluth_separation_{selected_target_shortname}.pdf"
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
            epsilon = float(row["bc_epsilon"])
            sigma_R = float(row["bc_sigma_R"])
            sigma_R_err = float(row["bc_sigma_R_err"])
            xbj = float(row["bc_xbj"])
            q2 = float(row["bc_q2"])

            if bin_num not in bin_data:
                bin_data[bin_num] = {"epsilon": [],
                                     "sigma_R": [],
                                     "sigma_R_err": [],
                                     "xbj": [],
                                     "q2": [],}
            bin_data[bin_num]["epsilon"].append(epsilon)
            bin_data[bin_num]["sigma_R"].append(sigma_R)
            bin_data[bin_num]["sigma_R_err"].append(sigma_R_err)
            bin_data[bin_num]["xbj"].append(xbj)
            bin_data[bin_num]["q2"].append(q2)
            
        except KeyError:
            continue
        except ValueError:
            continue

def linear_fit(x, intercept, slope):
    return intercept + slope * x

# for bin_num, data in bin_data.items():

#     eps = np.array(data["epsilon"])
#     sig = np.array(data["sigma_R"])
#     err = np.array(data["sigma_R_err"])
#     xbjavg = np.mean(data["xbj"])
#     n_points = len(data["q2"])
#     q2avg = np.mean(data["q2"])
#     q2max = np.max(data["q2"])
#     q2min = np.min(data["q2"])

#     if n_points > 2:
#         q2_err = np.std(data["q2"], ddof = 1)
#     else:
#         q2_err = (q2max - q2min) / 2

#     popt, pcov = curve_fit(linear_fit, eps, sig, sigma=err, absolute_sigma=True)

#     intercept, slope = popt
#     int_err = np.sqrt(pcov[0,0])
#     slope_err = np.sqrt(pcov[1,1])

#     R = float(slope) / float(intercept)

#     R_err = np.sqrt((pcov[1,1] / intercept**2) + (slope**2 * pcov[0,0] / intercept**4) - (2 * slope * pcov[0,1] / intercept**3))

#     print(f"bin {bin_num}: xbj = {xbj:.3f}, q2 = {q2avg:.3f} ± {q2_err:.3f}\nσ_T = {intercept:.4f} ± {int_err:.4f}, σ_L = {slope:.4f} ± {slope_err:.4f}, R = {R:.4f} ± {R_err:.4f}")
            
# plots_per_page = 2
# nrows, ncols = 2, 1
# fig = None
# plot_count = 0

# for bin_num, data in bin_data.items():
#     eps = np.array(data["epsilon"])
#     sig = np.array(data["sigma_R"])
#     err = np.array(data["sigma_R_err"])
#     xbjavg = np.mean(data["xbj"])
#     n_points = len(data["q2"])
#     q2avg = np.mean(data["q2"])
#     q2max = np.max(data["q2"])
#     q2min = np.min(data["q2"])
#     q2_err = np.std(data["q2"], ddof=1) if n_points > 2 else (q2max - q2min)/2

#     # linear fit
#     popt, pcov = curve_fit(linear_fit, eps, sig, sigma=err, absolute_sigma=True)
#     intercept, slope = popt

#     # create figure
#     fig, ax = plt.subplots(figsize=(6, 5))
    
#     # plot data + fit
#     ax.errorbar(eps, sig, yerr=err, fmt='o', color='black', label='data')
#     eps_fit = np.linspace(min(eps), max(eps), 100)
#     sig_fit = linear_fit(eps_fit, intercept, slope)
#     ax.plot(eps_fit, sig_fit, '-', color='red', label='fit')

#     # titles, labels
#     ax.set_title(f'Bin {bin_num} — xbj={xbjavg:.3f}, Q²={q2avg:.3f} ± {q2_err:.3f}', fontsize=10)
#     ax.set_xlabel('ε')
#     ax.set_ylabel('σ_R')
#     ax.legend(fontsize=8)
#     ax.grid(True)

#     # save page
#     pp.savefig(fig)
#     plt.close(fig)

# pp.close()
# print(f"PDF saved to {pdf_output}")

fit_results = []

with PdfPages(pdf_output) as pp:
    print(f"Found {len(bin_data)} bins. Processing...")

    for bin_num, data in bin_data.items():
        eps = np.array(data["epsilon"])
        sigR = np.array(data["sigma_R"])
        sigR_err = np.array(data["sigma_R_err"])

        eps_unique = {}
        for e, s, serr in zip(eps, sigR, sigR_err):
            if e not in eps_unique:
                eps_unique[e] = {"sigR": [], "sigR_err": []}
            eps_unique[e]["sigR"].append(s)
            eps_unique[e]["sigR_err"].append(serr)

        eps = []
        sigR = []
        sigR_err = []

        for e, vals in eps_unique.items():
            s = np.array(vals["sigR"])
            serr = np.array(vals["sigR_err"])

            weights = 1 / serr**2
            s_avg = np.sum(s * weights) / np.sum(weights)
            s_err = np.sqrt(1 / np.sum(weights))

            eps.append(e)
            sigR.append(s_avg)
            sigR_err.append(s_err)

        eps = np.array(eps)
        sigR = np.array(sigR)
        sigR_err = np.array(sigR_err)

        xbjavg = np.mean(data["xbj"])
        n_points = len(data["q2"])
        q2avg = np.mean(data["q2"])
        q2max = np.max(data["q2"])
        q2min = np.min(data["q2"])
        q2_err = np.std(data["q2"], ddof=1) if n_points > 2 else (q2max - q2min)/2

        popt, pcov = curve_fit(linear_fit, eps, sigR, sigma=sigR_err, absolute_sigma=True)
        intercept, slope = popt
        int_err = np.sqrt(pcov[0,0])
        slope_err = np.sqrt(pcov[1,1])
        R = float(slope) / float(intercept)
        R_err = np.sqrt((pcov[1,1] / intercept**2) + (slope**2 * pcov[0,0] / intercept**4) - (2 * slope * pcov[0,1] / intercept**3))

        print(f"bin {bin_num}: xbj = {xbjavg:.3f}, q2 = {q2avg:.3f} ± {q2_err:.3f}")
        print(f"σ_T = {intercept:.4f} ± {int_err:.4f}, σ_L = {slope:.4f} ± {slope_err:.4f}, R = {R:.4f} ± {R_err:.4f}")

        fit_results.append({"bin_num": bin_num,
                            "xbj": xbjavg,
                            "q2_avg": q2avg,
                            "q2_err": q2_err,
                            "sigma_T": intercept,
                            "sigma_T_err": int_err,
                            "sigma_L": slope,
                            "sigma_L_err": slope_err,
                            "R": R,
                            "R_err": R_err,})
        
        fig, ax = plt.subplots(figsize=(6, 5))
        
        ax.errorbar(eps, sigR, yerr=sigR_err, fmt="o", color="navy", markersize = 4)

        # Fit evaluated across full x-axis [0, 1]
        x_full = np.linspace(0.0, 1.0, 200)
        y_full = linear_fit(x_full, intercept, slope)
        ax.plot(x_full,y_full,"--",color="red",label=("Linear fit: y = mx + b\n"
                                                      rf"$\sigma_L = {slope:.4f} \pm {slope_err:.4f}$" "\n"
                                                      rf"$\sigma_T = {intercept:.4f} \pm {int_err:.4f}$"),)

        ax.set_title(f"{selected_target_titlename} Rosenbluth Separation\n"r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f} ± {q2_err:.3f}\n R = {R:.4f} ± {R_err:.4f}", fontsize=10)
        ax.set_xlabel(r"$\epsilon$")
        ax.set_ylabel(r"$d \sigma / d \Omega / dE' / \Gamma$ ($\mu$b/sr)")
        ax.set_xlim(0, 1)
        ax.legend(fontsize=8)
        ax.grid(True)

        pp.savefig(fig)
        plt.close(fig)

print(f"PDF saved to {pdf_output}")
csv_output = f"CSVs/{selected_run_type.upper()}_rosenbluth_fit_{selected_target_shortname}.csv"
pd.DataFrame(fit_results).to_csv(csv_output, index=False)
print(f"CSV of fits saved to {csv_output}")
