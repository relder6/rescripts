#!/usr/bin/env python3

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)
from INIT.config import get_common_run_inputs

# -----------------------------------------------------
# Inputs, files
# -----------------------------------------------------
selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, \
    selected_target_titlename, selected_target_A, selected_target_Z = get_common_run_inputs()

os.makedirs("CSVs", exist_ok=True)
os.makedirs("PDFs", exist_ok=True)

bc_csv = f"CSVs/{selected_run_type.upper()}_bin_centered_{selected_target_shortname}.csv"
pdf_output = f"PDFs/{selected_run_type.upper()}_rosenbluth_separation_{selected_target_shortname}.pdf"

# -----------------------------------------------------
# Read CSV
# -----------------------------------------------------
df = pd.read_csv(bc_csv)

# Sanity: require needed columns
required_cols = [
    "bin_num",
    "epsilon", "sigma_R_orig", "sigma_R_orig_err",  # pre-bin-centering
    "bc_epsilon", "bc_sigma_R", "bc_sigma_R_err",   # bin-centered
    "xbj", "bc_xbj", "q2", "bc_q2"
]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in {bc_csv}: {missing}")

# -----------------------------------------------------
# Group data per bin
# -----------------------------------------------------
bin_data_bc = {}   # bin-centered
bin_data_raw = {}  # pre-bin-centering

for bin_num, sub in df.groupby("bin_num"):
    bin_data_bc[bin_num] = {
        "epsilon": sub["bc_epsilon"].to_numpy(),
        "sigma_R": sub["bc_sigma_R"].to_numpy(),
        "sigma_R_err": sub["bc_sigma_R_err"].to_numpy(),
        "xbj": sub["bc_xbj"].to_numpy(),
        "q2": sub["bc_q2"].to_numpy(),
    }
    bin_data_raw[bin_num] = {
        "epsilon": sub["epsilon"].to_numpy(),
        "sigma_R": sub["sigma_R_orig"].to_numpy(),
        "sigma_R_err": sub["sigma_R_orig_err"].to_numpy(),
        "xbj": sub["xbj"].to_numpy(),
        "q2": sub["q2"].to_numpy(),
    }

# -----------------------------------------------------
# Fit function
# -----------------------------------------------------
def linear_fit(eps, intercept, slope):
    return intercept + slope * eps

# -----------------------------------------------------
# Plotting
# -----------------------------------------------------
from collections import defaultdict

with PdfPages(pdf_output) as pp:
    for bin_num in sorted(bin_data_bc.keys()):
        # ----- bin-centered per-event -----
        eps_bc = np.array(bin_data_bc[bin_num]["epsilon"])
        sig_bc = np.array(bin_data_bc[bin_num]["sigma_R"])
        err_bc = np.array(bin_data_bc[bin_num]["sigma_R_err"])

        # ----- pre-bin-centering per-event -----
        eps_raw = np.array(bin_data_raw[bin_num]["epsilon"])
        sig_raw = np.array(bin_data_raw[bin_num]["sigma_R"])
        err_raw = np.array(bin_data_raw[bin_num]["sigma_R_err"])

        # kinematics summary (use bin-centered q2/xbj)
        xbj_arr = np.array(bin_data_bc[bin_num]["xbj"])
        q2_arr = np.array(bin_data_bc[bin_num]["q2"])
        xbjavg = np.mean(xbj_arr)
        n_points = len(q2_arr)
        q2avg = np.mean(q2_arr)
        q2max = np.max(q2_arr)
        q2min = np.min(q2_arr)
        q2_err = np.std(q2_arr, ddof=1) if n_points > 2 else (q2max - q2min) / 2

        # ----- weighted averages per ε (bin-centered) -----
        eps_unique = defaultdict(lambda: {"sigR": [], "sigR_err": []})
        for e, s, serr in zip(eps_bc, sig_bc, err_bc):
            eps_unique[e]["sigR"].append(s)
            eps_unique[e]["sigR_err"].append(serr)

        eps_comb, sig_comb, err_comb = [], [], []
        for e, vals in eps_unique.items():
            s_arr = np.array(vals["sigR"])
            serr_arr = np.array(vals["sigR_err"])
            w = 1.0 / serr_arr**2
            s_avg = np.sum(s_arr * w) / np.sum(w)
            s_err = np.sqrt(1.0 / np.sum(w))
            eps_comb.append(e)
            sig_comb.append(s_avg)
            err_comb.append(s_err)

        eps_comb = np.array(eps_comb)
        sig_comb = np.array(sig_comb)
        err_comb = np.array(err_comb)

        # ----- fit to weighted bin-centered points -----
        popt, pcov = curve_fit(
            linear_fit, eps_comb, sig_comb, sigma=err_comb, absolute_sigma=True
        )
        intercept, slope = popt
        int_err = np.sqrt(pcov[0, 0])
        slope_err = np.sqrt(pcov[1, 1])
        R = slope / intercept
        R_err = np.sqrt(
            (pcov[1, 1] / intercept**2)
            + (slope**2 * pcov[0, 0] / intercept**4)
            - (2 * slope * pcov[0, 1] / intercept**3)
        )

        # ----- plot everything on one pad -----
        fig, ax = plt.subplots(figsize=(6, 5))

        # pre-bin-centering cloud
        ax.errorbar(
            eps_raw, sig_raw, yerr=err_raw,
            fmt="o", color="gray", alpha=0.5, markersize=3,
            label="pre bin-centering"
        )

        # bin-centered per-event points
        ax.errorbar(
            eps_bc, sig_bc, yerr=err_bc,
            fmt="o", color="lightgray", markersize=3,
            label="bin-centered points"
        )

        # weighted averages
        ax.errorbar(
            eps_comb, sig_comb, yerr=err_comb,
            fmt="o", color="navy", markersize=5,
            label="weighted avg per ε"
        )

        # fit line
        eps_min = min(eps_raw.min(), eps_bc.min(), eps_comb.min())
        eps_max = max(eps_raw.max(), eps_bc.max(), eps_comb.max())
        x_fit = np.linspace(eps_min, eps_max, 200)
        y_fit = linear_fit(x_fit, intercept, slope)
        ax.plot(
            x_fit, y_fit, "--", color="red",
            label=("Linear fit\n"
                   rf"$\sigma_L = {slope:.4f} \pm {slope_err:.4f}$" "\n"
                   rf"$\sigma_T = {intercept:.4f} \pm {int_err:.4f}$" "\n"
                   rf"$R = {R:.4f} \pm {R_err:.4f}$")
        )

        # axes limits
        all_sig = np.concatenate([sig_raw, sig_bc, sig_comb])
        y_min, y_max = all_sig.min(), all_sig.max()
        margin = 0.1
        dy = margin * (y_max - y_min if y_max > y_min else 1.0)
        ax.set_xlim(eps_min - 0.05, eps_max + 0.05)
        ax.set_ylim(y_min - dy, y_max + dy)

        # labels, title
        ax.set_title(
            f"{selected_target_titlename} Rosenbluth Separation\n"
            r"x$_{bj}$="f"{xbjavg:.3f}, Q"r"$^2$"f"={q2avg:.3f} ± {q2_err:.3f}",
            fontsize=10
        )
        ax.set_xlabel(r"$\epsilon$")
        ax.set_ylabel(r"$d \sigma / d \Omega / dE' / \Gamma$ ($\mu$b/sr)")
        ax.legend(fontsize=8)
        ax.grid(True)

        pp.savefig(fig)
        plt.close(fig)

print(f"PDF saved to {pdf_output}")
