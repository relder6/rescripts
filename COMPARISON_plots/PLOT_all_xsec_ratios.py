#!/usr/bin/env python3

import os, sys, csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
import mplhep
from collections import defaultdict
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import parse_beam_pass, get_common_values

csv_file = "CSVs/RME_xsec_ratio_results.csv"

arg1 = sys.argv[1] if len(sys.argv) > 1 else None
selected_beam_pass, beam_prefix = parse_beam_pass(arg1)

vals = get_common_values()

if selected_beam_pass == "4":
    ebeam = float(vals["ebeam_4pass"])
    theta = float(vals["angle_4pass"])
elif selected_beam_pass == "5":
    ebeam = float(vals["ebeam_5pass"])
    theta = float(vals["angle_5pass"])

if not os.path.exists(csv_file):
    print(f"Missing input file: {csv_file}.  Try again.")
    exit

data = np.genfromtxt(csv_file, delimiter = ",", names=True, dtype=None, encoding=None)

pass_mask = np.isin(data["exp"], [f"RSIDIS({selected_beam_pass}Pass)"])

targs = ["Aluminum", "Carbon", "Copper"]

# targs = np.unique(data["target"][pass_mask])

plt.style.use(mplhep.style.ROOT)

plt.rcParams.update({"font.family": "DejaVu Sans",
                     "mathtext.fontset": "dejavusans",
                     "mathtext.default": "regular",
                     "figure.titlesize": 14,
                     "axes.titlesize": 16,
                     "axes.labelsize": 14,
                     "legend.fontsize": 12,
                     "xtick.labelsize": 12,
                     "ytick.labelsize": 12,})

        
fig, ax = plt.subplots(figsize=(6.5, 5),constrained_layout=True)
for targ in targs:
    targ_mask = data["target_num"] == targ
    mask = pass_mask & targ_mask
    xsec = data["xsec_ratio_per_nucleon"][mask]
    xsec_err = data["xsec_ratio_per_nucleon_err"][mask]
    eprime = data["eprime"][mask]
    q2 = data["q2"][mask]
    xbj = data["xbj"][mask]
    nu = (ebeam - data["eprime"])[mask]
    print("selected_beam_pass =", selected_beam_pass)
    print("unique exp values =", np.unique(data["exp"]))
    print("matched rows =", np.sum(pass_mask))

    ax.errorbar(q2, xsec, yerr=xsec_err, fmt = "o", markersize = 4, capsize = 0, label = f"{targ}")
    
ax.set_ylabel(r"$\sigma_A / A / \sigma_D / D$")
ax.set_title(f"Cross Section Ratios per Nucleon\nE = {ebeam:.2f} GeV at {theta:.2f}$^\circ$")
ax.set_xlabel(r"$Q^2$ [$GeV^2$]", labelpad=-2, loc="right")

q2_ticks = np.linspace(*ax.get_xlim(), 6)

ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
ax.set_xticks(q2_ticks)

sort_idx = np.argsort(q2)
nu_ticks = np.interp(q2_ticks, q2[sort_idx], nu[sort_idx])

ax_bot = ax.secondary_xaxis(-0.25)
ax_bot.set_xlabel(r"$\nu$ [$GeV$]", labelpad=-2, loc="right")
ax_bot.set_xticks(q2_ticks)
ax_bot.set_xticklabels([f"{x:.2f}" for x in nu_ticks])

ax_bot.tick_params(axis="x", which="both", direction="in", pad=5, length=4)
ax_bot.tick_params(axis="x", which="major", direction="in", pad=5, length=10, width=1)

ax.grid()
ax.legend(loc="upper right",frameon=True,fancybox=True,framealpha=0.6,edgecolor="gray")
plt.show()
