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
from INIT.config import get_common_run_inputs
from scipy.optimize import curve_fit
from collections import defaultdict

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, selected_target_titlename, selected_target_A, selected_target_Z = get_common_run_inputs()

os.makedirs("CSVs", exist_ok=True)
bc_csv = f"CSVs/{selected_run_type.upper()}_bin_centered_{selected_target_shortname}.csv"

os.makedirs("PDFs", exist_ok=True)
pdf_output = f"PDFs/{selected_run_type.upper()}_rosenbluth_separation_{selected_target_shortname}.pdf"
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

# fit_results = []
# eps_rows = []

# with PdfPages(pdf_output) as pp:
#     print(f"Found {len(bin_data)} bins. Processing...")

#     for bin_num, data in bin_data.items():
#         eps = np.array(data["epsilon"])
#         sigR = np.array(data["sigma_R"])
#         sigR_err = np.array(data["sigma_R_err"])

#         eps_unique = {}
#         for e, s, serr in zip(eps, sigR, sigR_err):
#             if e not in eps_unique:
#                 eps_unique[e] = {"sigR": [], "sigR_err": []}
#             eps_unique[e]["sigR"].append(s)
#             eps_unique[e]["sigR_err"].append(serr)

#         eps = []
#         sigR = []
#         sigR_err = []

#         for e, vals in eps_unique.items():
#             s = np.array(vals["sigR"])
#             serr = np.array(vals["sigR_err"])

#             weights = 1 / serr**2
#             s_avg = np.sum(s * weights) / np.sum(weights)
#             s_err = np.sqrt(1 / np.sum(weights))

#             eps.append(e)
#             sigR.append(s_avg)
#             sigR_err.append(s_err)

#         eps = np.array(eps)
#         sigR = np.array(sigR)
#         sigR_err = np.array(sigR_err)

#         for e, s, serr in zip(eps, sigR, sigR_err):
#             eps_rows.append({"bin_num": bin_num,
#                              "epsilon": e,
#                              "sigma_R": s,
#                              "sigma_R_err": serr})

#         xbjavg = np.mean(data["xbj"])
#         n_points = len(data["q2"])
#         q2avg = np.mean(data["q2"])
#         q2max = np.max(data["q2"])
#         q2min = np.min(data["q2"])
#         q2_err = np.std(data["q2"], ddof=1) if n_points > 2 else (q2max - q2min)/2

#         print(f"\nBIN {bin_num}")
#         print("eps:", eps)
#         print("sigR:", sigR)
#         print("sigR_err:", sigR_err)
#         print("N points:", len(eps))

#         popt, pcov = curve_fit(linear_fit, eps, sigR, sigma=sigR_err, absolute_sigma=True)
#         intercept, slope = popt
#         int_err = np.sqrt(pcov[0,0])
#         slope_err = np.sqrt(pcov[1,1])
#         R = float(slope) / float(intercept)
#         R_err = np.sqrt((pcov[1,1] / intercept**2) + (slope**2 * pcov[0,0] / intercept**4) - (2 * slope * pcov[0,1] / intercept**3))

#         print(f"bin {bin_num}: xbj = {xbjavg:.3f}, q2 = {q2avg:.3f} ± {q2_err:.3f}")
#         print(f"σ_T = {intercept:.4f} ± {int_err:.4f}, σ_L = {slope:.4f} ± {slope_err:.4f}, R = {R:.4f} ± {R_err:.4f}")

#         for e, s, serr in zip(eps, sigR, sigR_err):
#             fit_results.append({"bin_num": bin_num,
#                                 "xbj": xbjavg,
#                                 "q2_avg": q2avg,
#                                 "q2_err": q2_err,
#                                 "epsilon": e,
#                                 "sigma_R": s,
#                                 "sigma_R_err": serr,
#                                 "sigma_T": intercept,
#                                 "sigma_T_err": int_err,
#                                 "sigma_L": slope,
#                                 "sigma_L_err": slope_err,
#                                 "R": R,
#                                 "R_err": R_err})
        
#         fig, ax = plt.subplots(figsize=(6, 5))
        
#         ax.errorbar(eps, sigR, yerr=sigR_err, fmt="o", color="navy", markersize = 4, label = "test")

#         margin = 0.1 #5% plotting margin
        
#         x_margin = margin * (max(eps) - min(eps))
#         y_margin = margin * (max(sigR) - min(sigR))

#         ax.set_xlim(min(eps) - x_margin, max(eps) + x_margin)
#         ax.set_ylim(min(sigR) - y_margin, max(sigR) + y_margin)
        
#         x_full = np.linspace(min(eps) - x_margin, max(eps) + x_margin, 200)
        
#         y_full = linear_fit(x_full, intercept, slope)
#         ax.plot(x_full,y_full,"--",color="red",label=("Linear fit: y = mx + b\n"
#                                                       rf"$\sigma_L = {slope:.4f} \pm {slope_err:.4f}$" "\n"
#                                                       rf"$\sigma_T = {intercept:.4f} \pm {int_err:.4f}$"),)

#         ax.set_title(f"{selected_target_titlename} Rosenbluth Separation\n"r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f} ± {q2_err:.3f}\n R = {R:.4f} ± {R_err:.4f}", fontsize=10)
#         ax.set_xlabel(r"$\epsilon$")
#         ax.set_ylabel(r"$d \sigma / d \Omega / dE' / \Gamma$ ($\mu$b/sr)")
#         ax.legend(fontsize=8)
#         ax.grid(True)

#         pp.savefig(fig)
#         plt.close(fig)

# with PdfPages(pdf_output) as pp:
#     print(f"Found {len(bin_data)} bins. Processing...")

#     for bin_num, data in bin_data.items():
#         # --------------------------
#         # Raw points
#         # --------------------------
#         eps_all = np.array(data["epsilon"])
#         sigR_all = np.array(data["sigma_R"])
#         sigR_err_all = np.array(data["sigma_R_err"])

#         # --------------------------
#         # Weighted-average points
#         # --------------------------
#         eps_unique = {}
#         for e, s, serr in zip(eps_all, sigR_all, sigR_err_all):
#             if e not in eps_unique:
#                 eps_unique[e] = {"sigR": [], "sigR_err": []}
#             eps_unique[e]["sigR"].append(s)
#             eps_unique[e]["sigR_err"].append(serr)

#         eps_w = []
#         sigR_w = []
#         sigR_err_w = []

#         for e, vals in eps_unique.items():
#             s = np.array(vals["sigR"])
#             serr = np.array(vals["sigR_err"])
#             weights = 1 / serr**2
#             s_avg = np.sum(s * weights) / np.sum(weights)
#             s_err = np.sqrt(1 / np.sum(weights))
#             eps_w.append(e)
#             sigR_w.append(s_avg)
#             sigR_err_w.append(s_err)

#         eps_w = np.array(eps_w)
#         sigR_w = np.array(sigR_w)
#         sigR_err_w = np.array(sigR_err_w)

#         # --------------------------
#         # Fit to weighted points
#         # --------------------------
#         popt_w, pcov_w = curve_fit(linear_fit, eps_w, sigR_w, sigma=sigR_err_w, absolute_sigma=True)
#         intercept_w, slope_w = popt_w
#         int_err_w = np.sqrt(pcov_w[0,0])
#         slope_err_w = np.sqrt(pcov_w[1,1])

#         # --------------------------
#         # Fit to all points
#         # --------------------------
#         popt_all, pcov_all = curve_fit(linear_fit, eps_all, sigR_all, sigma=sigR_err_all, absolute_sigma=True)
#         intercept_all, slope_all = popt_all
#         int_err_all = np.sqrt(pcov_all[0,0])
#         slope_err_all = np.sqrt(pcov_all[1,1])

#         # --------------------------
#         # Plot 1: weighted points
#         # --------------------------
#         fig, ax = plt.subplots(figsize=(6,5))
#         ax.errorbar(eps_w, sigR_w, yerr=sigR_err_w, fmt='o', color='navy', label='Weighted points')
#         x_fit = np.linspace(min(eps_w), max(eps_w), 200)
#         ax.plot(x_fit, linear_fit(x_fit, intercept_w, slope_w), '--', color='red',
#                 label=(f'Weighted Fit:\nσ_L = {slope_w:.4f} ± {slope_err_w:.4f}\nσ_T = {intercept_w:.4f} ± {int_err_w:.4f}'))
#         ax.set_title(f"Bin {bin_num} — Weighted points")
#         ax.set_xlabel("ε")
#         ax.set_ylabel(r"$\sigma_R$")
#         ax.grid(True)
#         ax.legend(fontsize=8)
#         pp.savefig(fig)
#         plt.close(fig)

#         # --------------------------
#         # Plot 2: all raw points
#         # --------------------------
#         fig, ax = plt.subplots(figsize=(6,5))
#         ax.errorbar(eps_all, sigR_all, yerr=sigR_err_all, fmt='o', color='darkgreen', label='All points')
#         x_fit = np.linspace(min(eps_all), max(eps_all), 200)
#         ax.plot(x_fit, linear_fit(x_fit, intercept_all, slope_all), '--', color='red',
#                 label=(f'Fit to all points:\nσ_L = {slope_all:.4f} ± {slope_err_all:.4f}\nσ_T = {intercept_all:.4f} ± {int_err_all:.4f}'))
#         ax.set_title(f"Bin {bin_num} — All raw points")
#         ax.set_xlabel("ε")
#         ax.set_ylabel(r"$\sigma_R$")
#         ax.grid(True)
#         ax.legend(fontsize=8)
#         pp.savefig(fig)
#         plt.close(fig)

# ~~~~~THis one is pretty good

# from collections import defaultdict
# from matplotlib.backends.backend_pdf import PdfPages
# from scipy.optimize import curve_fit

# def linear_fit(eps, intercept, slope):
#     return intercept + slope * eps

# with PdfPages(pdf_output) as pp:
#     print(f"Found {len(bin_data)} bins. Processing...")

#     for bin_num, data in bin_data.items():
#         # ---------- Common bin metadata ----------
#         eps_raw = np.array(data["epsilon"])
#         sig_raw = np.array(data["sigma_R"])
#         err_raw = np.array(data["sigma_R_err"])

#         xbjavg = np.mean(data["xbj"])
#         n_points = len(data["q2"])
#         q2avg = np.mean(data["q2"])
#         q2max = np.max(data["q2"])
#         q2min = np.min(data["q2"])
#         q2_err = np.std(data["q2"], ddof=1) if n_points > 2 else (q2max - q2min) / 2

#         # ---------- RAW FIT ----------
#         popt_raw, pcov_raw = curve_fit(linear_fit, eps_raw, sig_raw,
#                                        sigma=err_raw, absolute_sigma=True)
#         intercept_raw, slope_raw = popt_raw
#         int_err_raw = np.sqrt(pcov_raw[0, 0])
#         slope_err_raw = np.sqrt(pcov_raw[1, 1])
#         R_raw = slope_raw / intercept_raw
#         R_err_raw = np.sqrt(
#             (pcov_raw[1, 1] / intercept_raw**2) +
#             (slope_raw**2 * pcov_raw[0, 0] / intercept_raw**4) -
#             (2 * slope_raw * pcov_raw[0, 1] / intercept_raw**3)
#         )

#         # ---------- WEIGHTED-MEAN BY EPSILON ----------
#         eps_unique = defaultdict(lambda: {"sigR": [], "sigR_err": []})
#         for e, s, serr in zip(eps_raw, sig_raw, err_raw):
#             eps_unique[e]["sigR"].append(s)
#             eps_unique[e]["sigR_err"].append(serr)

#         eps_comb, sig_comb, err_comb = [], [], []
#         for e, vals in eps_unique.items():
#             s_arr = np.array(vals["sigR"])
#             serr_arr = np.array(vals["sigR_err"])
#             weights = 1.0 / serr_arr**2
#             s_avg = np.sum(s_arr * weights) / np.sum(weights)
#             s_err = np.sqrt(1.0 / np.sum(weights))
#             eps_comb.append(e)
#             sig_comb.append(s_avg)
#             err_comb.append(s_err)

#         eps_comb = np.array(eps_comb)
#         sig_comb = np.array(sig_comb)
#         err_comb = np.array(err_comb)

#         popt_comb, pcov_comb = curve_fit(linear_fit, eps_comb, sig_comb,
#                                          sigma=err_comb, absolute_sigma=True)
#         intercept_comb, slope_comb = popt_comb
#         int_err_comb = np.sqrt(pcov_comb[0, 0])
#         slope_err_comb = np.sqrt(pcov_comb[1, 1])
#         R_comb = slope_comb / intercept_comb
#         R_err_comb = np.sqrt(
#             (pcov_comb[1, 1] / intercept_comb**2) +
#             (slope_comb**2 * pcov_comb[0, 0] / intercept_comb**4) -
#             (2 * slope_comb * pcov_comb[0, 1] / intercept_comb**3)
#         )

#         # ---------- PLOTTING SIDE BY SIDE ----------
#         fig, (ax_raw, ax_comb) = plt.subplots(1, 2, figsize=(11, 4), sharey=True)

#         # margins for consistent axes
#         def add_margins(x, y, margin=0.1):
#             x_min, x_max = x.min(), x.max()
#             y_min, y_max = y.min(), y.max()
#             dx = margin * (x_max - x_min if x_max > x_min else 1.0)
#             dy = margin * (y_max - y_min if y_max > y_min else 1.0)
#             return (x_min - dx, x_max + dx, y_min - dy, y_max + dy)

#         # ----- RAW PANEL -----
#         ax_raw.errorbar(eps_raw, sig_raw, yerr=err_raw, fmt="o",
#                         color="black", markersize=4, label="raw data")
#         x1, x2, y1, y2 = add_margins(eps_raw, sig_raw)
#         ax_raw.set_xlim(x1, x2)
#         ax_raw.set_ylim(y1, y2)

#         x_full_raw = np.linspace(x1, x2, 200)
#         y_full_raw = linear_fit(x_full_raw, intercept_raw, slope_raw)
#         ax_raw.plot(
#             x_full_raw, y_full_raw, "--", color="red",
#             label=("Linear fit\n"
#                    rf"$\sigma_L = {slope_raw:.4f} \pm {slope_err_raw:.4f}$" "\n"
#                    rf"$\sigma_T = {intercept_raw:.4f} \pm {int_err_raw:.4f}$" "\n"
#                    rf"$R = {R_raw:.4f} \pm {R_err_raw:.4f}$")
#         )

#         ax_raw.set_title(
#             f"{selected_target_titlename} (raw)\n"
#             r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f} ± {q2_err:.3f}",
#             fontsize=9
#         )
#         ax_raw.set_xlabel(r"$\epsilon$")
#         ax_raw.set_ylabel(r"$d \sigma / d \Omega / dE' / \Gamma$ ($\mu$b/sr)")
#         ax_raw.legend(fontsize=7)
#         ax_raw.grid(True)

#         print(rf"$R = {R_raw:.4f} \pm {R_err_raw:.4f}$")

#         # ----- COMBINED PANEL -----
#         ax_comb.errorbar(eps_comb, sig_comb, yerr=err_comb, fmt="o",
#                          color="navy", markersize=4, label="weighted means")
#         x1c, x2c, y1c, y2c = add_margins(eps_comb, sig_comb)
#         ax_comb.set_xlim(x1c, x2c)
#         ax_comb.set_ylim(y1, y2)  # share y-range with raw for visual comparison

#         x_full_comb = np.linspace(x1c, x2c, 200)
#         y_full_comb = linear_fit(x_full_comb, intercept_comb, slope_comb)
#         ax_comb.plot(
#             x_full_comb, y_full_comb, "--", color="red",
#             label=("Linear fit\n"
#                    rf"$\sigma_L = {slope_comb:.4f} \pm {slope_err_comb:.4f}$" "\n"
#                    rf"$\sigma_T = {intercept_comb:.4f} \pm {int_err_comb:.4f}$" "\n"
#                    rf"$R = {R_comb:.4f} \pm {R_err_comb:.4f}$")
#         )

#         ax_comb.set_title(
#             f"{selected_target_titlename} (Weighted Avg)\n"
#             r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f} ± {q2_err:.3f}",
#             fontsize=9
#         )
#         ax_comb.set_xlabel(r"$\epsilon$")
#         ax_comb.legend(fontsize=7)
#         ax_comb.grid(True)

#         fig.suptitle(f"Bin {bin_num} Rosenbluth Separation", fontsize=10)
#         fig.tight_layout(rect=[0, 0.03, 1, 0.95])
#         pp.savefig(fig)
#         plt.close(fig)

# from collections import defaultdict
# from scipy.optimize import curve_fit

# def linear_fit(eps, intercept, slope):
#     return intercept + slope * eps

# with PdfPages(pdf_output) as pp:
#     for bin_num, data in bin_data.items():
#         # ----- RAW (bin-centered) points -----
#         eps_raw = np.array(data["epsilon"])       # or bc_epsilon
#         sig_raw = np.array(data["sigma_R"])       # or bc_sigma_R
#         err_raw = np.array(data["sigma_R_err"])   # or bc_sigma_R_err

#         xbjavg = np.mean(data["xbj"])
#         n_points = len(data["q2"])
#         q2avg = np.mean(data["q2"])
#         q2max = np.max(data["q2"])
#         q2min = np.min(data["q2"])
#         q2_err = np.std(data["q2"], ddof=1) if n_points > 2 else (q2max - q2min)/2

#         # ----- RAW FIT -----
#         popt_raw, pcov_raw = curve_fit(
#             linear_fit, eps_raw, sig_raw, sigma=err_raw, absolute_sigma=True
#         )
#         intercept_raw, slope_raw = popt_raw
#         int_err_raw = np.sqrt(pcov_raw[0, 0])
#         slope_err_raw = np.sqrt(pcov_raw[1, 1])
#         R_raw = slope_raw / intercept_raw
#         R_err_raw = np.sqrt(
#             (pcov_raw[1, 1] / intercept_raw**2) +
#             (slope_raw**2 * pcov_raw[0, 0] / intercept_raw**4) -
#             (2 * slope_raw * pcov_raw[0, 1] / intercept_raw**3)
#         )

#         # ----- BUILD WEIGHTED AVERAGES AT SAME ε -----
#         eps_unique = defaultdict(lambda: {"sigR": [], "sigR_err": []})
#         for e, s, serr in zip(eps_raw, sig_raw, err_raw):
#             eps_unique[e]["sigR"].append(s)
#             eps_unique[e]["sigR_err"].append(serr)

#         eps_comb, sig_comb, err_comb = [], [], []
#         for e, vals in eps_unique.items():
#             s_arr = np.array(vals["sigR"])
#             serr_arr = np.array(vals["sigR_err"])
#             w = 1.0 / serr_arr**2
#             s_avg = np.sum(s_arr * w) / np.sum(w)
#             s_err = np.sqrt(1.0 / np.sum(w))
#             eps_comb.append(e)
#             sig_comb.append(s_avg)
#             err_comb.append(s_err)

#         eps_comb = np.array(eps_comb)
#         sig_comb = np.array(sig_comb)
#         err_comb = np.array(err_comb)

#         # Fit on the *combined* points (will numerically match raw fit, but we want them drawn)
#         popt_comb, pcov_comb = curve_fit(
#             linear_fit, eps_comb, sig_comb, sigma=err_comb, absolute_sigma=True
#         )
#         intercept_comb, slope_comb = popt_comb
#         int_err_comb = np.sqrt(pcov_comb[0, 0])
#         slope_err_comb = np.sqrt(pcov_comb[1, 1])
#         R_comb = slope_comb / intercept_comb
#         R_err_comb = np.sqrt(
#             (pcov_comb[1, 1] / intercept_comb**2) +
#             (slope_comb**2 * pcov_comb[0, 0] / intercept_comb**4) -
#             (2 * slope_comb * pcov_comb[0, 1] / intercept_comb**3)
#         )

#         # ----- SINGLE PAD: raw + weighted-average points -----
#         fig, ax = plt.subplots(figsize=(6, 5))

#         # raw points
#         ax.errorbar(
#             eps_raw, sig_raw, yerr=err_raw,
#             fmt="o", color="lightgray", markersize=3,
#             label="raw bin-centered points"
#         )

#         # weighted-average points (emphasized)
#         ax.errorbar(
#             eps_comb, sig_comb, yerr=err_comb,
#             fmt="o", color="navy", markersize=5,
#             label="weighted avg per ε"
#         )

#         # Draw only one fit line (they’re the same numerically). Use combined params.
#         eps_min = min(eps_raw.min(), eps_comb.min())
#         eps_max = max(eps_raw.max(), eps_comb.max())
#         x_fit = np.linspace(eps_min, eps_max, 200)
#         y_fit = linear_fit(x_fit, intercept_comb, slope_comb)
#         ax.plot(
#             x_fit, y_fit, "--", color="red",
#             label=("Linear fit\n"
#                    rf"$\sigma_L = {slope_comb:.4f} \pm {slope_err_comb:.4f}$" "\n"
#                    rf"$\sigma_T = {intercept_comb:.4f} \pm {int_err_comb:.4f}$" "\n"
#                    rf"$R = {R_comb:.4f} \pm {R_err_comb:.4f}$")
#         )

#         # cosmetic stuff
#         margin = 0.1
#         y_min, y_max = sig_raw.min(), sig_raw.max()
#         dy = margin * (y_max - y_min if y_max > y_min else 1.0)
#         ax.set_xlim(eps_min - 0.05, eps_max + 0.05)
#         ax.set_ylim(y_min - dy, y_max + dy)

#         ax.set_title(
#             f"{selected_target_titlename} Rosenbluth Separation\n"
#             r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f} ± {q2_err:.3f}",
#             fontsize=10
#         )
#         ax.set_xlabel(r"$\epsilon$")
#         ax.set_ylabel(r"$d \sigma / d \Omega / dE' / \Gamma$ ($\mu$b/sr)")
#         ax.legend(fontsize=8)
#         ax.grid(True)

#         pp.savefig(fig)
#         plt.close(fig)

# print(f"PDF saved to {pdf_output}")

# # Example: building a second dict for "raw" epsilon & sigma_R
# bin_data_raw = {}
# for bin_num in sorted(df_final["bin_num"].unique()):
#     sub = df_final[df_final["bin_num"] == bin_num]
#     bin_data_raw[bin_num] = {
#         "epsilon": sub["epsilon"].to_numpy(),            # pre-bin-centering ε
#         "sigma_R": sub["sigma_R"].to_numpy(),            # pre-bin-centering σ_R
#         "sigma_R_err": sub["sigma_R_err"].to_numpy(),
#         "xbj": sub["xbj"].to_numpy(),
#         "q2": sub["q2"].to_numpy(),
#     }
# with PdfPages(pdf_smear_output) as pp:
#     for bin_num in bin_data.keys():
#         # after centering (current bin_data)
#         eps_bc  = np.array(bin_data[bin_num]["epsilon"])       # bc_epsilon
#         sig_bc  = np.array(bin_data[bin_num]["sigma_R"])       # bc_sigma_R
#         err_bc  = np.array(bin_data[bin_num]["sigma_R_err"])

#         # before centering
#         eps_raw = np.array(bin_data_raw[bin_num]["epsilon"])
#         sig_raw = np.array(bin_data_raw[bin_num]["sigma_R"])
#         err_raw = np.array(bin_data_raw[bin_num]["sigma_R_err"])

#         xbjavg = np.mean(bin_data_raw[bin_num]["xbj"])
#         q2avg  = np.mean(bin_data_raw[bin_num]["q2"])

#         fig, ax = plt.subplots(figsize=(6, 5))

#         # pre-bin-centering: show smear
#         ax.errorbar(
#             eps_raw, sig_raw, yerr=err_raw,
#             fmt="o", color="gray", alpha=0.6, markersize=3,
#             label="pre-bin-centering"
#         )

#         # post-bin-centering: show tightened points
#         ax.errorbar(
#             eps_bc, sig_bc, yerr=err_bc,
#             fmt="o", color="darkgreen", markersize=4,
#             label="post-bin-centering"
#         )

#         # optional: arrows from raw → centered for a few points (if 1–1 mapping exists)
#         # for e0, s0, e1, s1 in zip(eps_raw, sig_raw, eps_bc, sig_bc):
#         #     ax.arrow(e0, s0, e1 - e0, s1 - s0,
#         #              length_includes_head=True,
#         #              head_width=0.01, head_length=0.02,
#         #              color="k", alpha=0.3)

#         margin = 0.1
#         all_eps = np.concatenate([eps_raw, eps_bc])
#         all_sig = np.concatenate([sig_raw, sig_bc])
#         e_min, e_max = all_eps.min(), all_eps.max()
#         s_min, s_max = all_sig.min(), all_sig.max()
#         de = margin * (e_max - e_min if e_max > e_min else 1.0)
#         ds = margin * (s_max - s_min if s_max > s_min else 1.0)
#         ax.set_xlim(e_min - de, e_max + de)
#         ax.set_ylim(s_min - ds, s_max + ds)

#         ax.set_title(
#             f"Bin {bin_num}: pre vs post bin-centering\n"
#             r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f}",
#             fontsize=10
#         )
#         ax.set_xlabel(r"$\epsilon$")
#         ax.set_ylabel(r"$d \sigma / d \Omega / dE' / \Gamma$ ($\mu$b/sr)")
#         ax.legend(fontsize=8)
#         ax.grid(True)

#         pp.savefig(fig)
#         plt.close(fig)

with PdfPages(pdf_output) as pp:
    for bin_num in sorted(bc_data.keys()):

        eps_bc = np.array(bc_data[bin_num]["epsilon"])
        sig_bc = np.array(bc_data[bin_num]["sigma_R"])
        sig_err_bc = np.array(bc_data[bin_num]["sigma_R_err"])
        eps_raw = np.array(data[bin_num]["epsilon"])
        sig_raw = np.array(data[bin_num]["sigma_R"])
        sig_err_raw = np.array(data[bin_num]["sigma_R_err"])

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
                                                       rf"$R = {R:.4f} \pm {R_err:.4f}$"))

        all_sig = np.concatenate([sig_raw, sig_bc, sig_comb])
        y_min, y_max = all_sig.min(), all_sig.max()
        margin = 0.1
        dy = margin * (y_max - y_min if y_max > y_min else 1.0)
        ax.set_xlim(eps_min - 0.05, eps_max + 0.05)
        ax.set_ylim(y_min - dy, y_max + dy)

        ax.set_title(f"{selected_target_titlename} Rosenbluth Separation\n"
                     r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f} ± {q2_err:.3f}",fontsize=10)
        ax.set_xlabel(r"$\epsilon$")
        ax.set_ylabel(r"$d \sigma / d \Omega / dE' / \Gamma$ ($\mu$b/sr)")
        ax.legend(fontsize=8)
        ax.grid(True)

        pp.savefig(fig)
        plt.close(fig)

print(f"PDF saved to {pdf_output}")

csv_output = f"CSVs/{selected_run_type.upper()}_rosenbluth_fit_{selected_target_shortname}.csv"
pd.DataFrame(fit_results).to_csv(csv_output, index=False)
print(f"CSV of fits saved to {csv_output}")
