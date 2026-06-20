#!/usr/bin/env python3

import pandas as pd
import os, re, sys
import numpy as np
import matplotlib.ticker as ticker
import mplhep
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as mtick

# rel_guess = 0.75 / 100.0 # 0.5% from Deb's thesis, I'm assuming worse here.

info_directory = "../../FILTER_type"
output_csv = "CSVs/testing.csv"

targets = ["Al", "C", "Cu", "ld2", "lh2"]

beam_passes = ["4pass", "5pass"]

all_rows = []

for target in targets:
    for beam_pass in beam_passes:
        csv_file = f"{info_directory}/{target.upper()}/hmsdis_{beam_pass}_{target.lower()}_runs.csv"
        if not os.path.exists(csv_file):
            print(f"Skipping missing file: {csv_file}")
            continue

        df = pd.read_csv(csv_file)

        df["target"] = target
        df["beam_pass"] = beam_pass
        
        df_out = df[["target", "beam_pass", "runnum", "ibeam", "ibeam1", "ibeam4a", "ibeam4c"]]

        df_out = df_out.rename(columns={"ibeam": "ibeam2"})

        df_out = df_out.dropna()

        all_rows.append(df_out)
        
df_all = pd.concat(all_rows, ignore_index=True)

df_all["i2_ov_i1"] = df_all["ibeam2"] / df_all["ibeam1"]

df_all["i2_ov_i4a"] = df_all["ibeam2"] / df_all["ibeam4a"]

df_all["i2_ov_i4c"] = df_all["ibeam2"] / df_all["ibeam4c"]

df_all["i4a_ov_i4c"] = df_all["ibeam4a"] / df_all["ibeam4c"]

df_all["ibeam2_avg"] = (df_all.groupby(["target", "beam_pass"])["ibeam2"].transform("mean"))

df_all["syst_abs"] = 0.200 / df_all["ibeam2_avg"]

df_all["n_entries"] = (df_all.groupby(["target", "beam_pass"])["runnum"].transform("count"))

df_all["bcm2_inferred_rate"] = (df_all["ibeam2"] * 5677.2 + 250068)

df_all["bcm2_calib_1_i"] = (df_all["bcm2_inferred_rate"] - 250382) / 5677.1

df_all["bcm2_calib_2_i"] = (df_all["bcm2_inferred_rate"] - 250718) / 5645.7

df_all["bcm2_calib_i1_ov_i2"] = 1 - (df_all["bcm2_calib_1_i"] / df_all["bcm2_calib_2_i"])

# df_all["syst_rel"] = rel_guess / np.sqrt(df_all["n_entries"])

# df_all["syst_tot"] = np.sqrt((df_all["syst_rel"]**2) + (df_all["syst_abs"]**2))

df_all.to_csv(output_csv, index=False)

print(f"Saved compiled CSV → {output_csv}")

import matplotlib.pyplot as plt

# ----------------------------------
# Making a plot
# ----------------------------------
# df_plot = df_all.drop_duplicates(subset=["target", "beam_pass"])

# df_plot = df_plot.sort_values(["beam_pass", "target"])

# labels = df_plot["target"] + "_" + df_plot["beam_pass"]

# x = np.arange(len(labels))

# width = 0.25

# ----------------------------------
# Plots
# ----------------------------------
# plt.style.use(mplhep.style.ROOT)
# plt.rcParams.update({"figure.titlesize": 14,
#                      "axes.titlesize": 12,
#                      "axes.labelsize": 10,
#                      "legend.fontsize": 10,
#                      "xtick.labelsize": 10,
#                      "ytick.labelsize": 10})

# fig, ax = plt.subplots(figsize=(8,6))

# ax.bar(x - width, df_plot["syst_rel"]*100, width, label=f"syst_rel (assuming bcms tracking to {100*rel_guess:.2f}%)")
# ax.bar(x, df_plot["syst_abs"]*100, width, label="syst_abs")
# ax.bar(x + width, df_plot["syst_tot"]*100, width, label="syst_tot")

# ax.set_xticks(x)
# ax.set_xticklabels(labels, rotation=45)
# ax.set_ylabel("Systematic Uncertainty (Percent)")
# ax.set_title("Beam Charge Systematics by Target and Beam Pass")
# ax.yaxis.set_major_formatter(ticker.PercentFormatter(100))

# ax.legend()

# plt.grid(axis = "y")

# plt.tight_layout()
# plt.savefig("PNGs/beam_charge_systematics.png", dpi=300)
# plt.show()

# plt.close()

# ------------------------------------------------------
# Defining some things for plotting
# ------------------------------------------------------
# ibeams = ["ibeam1", "ibeam2", "ibeam4a", "ibeam4c"]
# ratios = ["i2_ov_i1", "i2_ov_i4a", "i2_ov_i4c", "i4a_ov_i4c"]

ibeams = ["ibeam2", "ibeam4a", "ibeam4c"]
ratios = ["i2_ov_i4a", "i2_ov_i4c", "i4a_ov_i4c"]
calibs = ["bcm2_calib_1_i", "bcm2_calib_2_i"]
beam_passes = ["4pass", "5pass"]

def p0_fit(y):
    y = np.asarray(y)
    p0 = np.mean(y)
    return p0

# ------------------------------------------------------
# First plot, bcm current and ratios vs run number
# ------------------------------------------------------
fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize = (8, 6), gridspec_kw = {"height_ratios":[1,4], "hspace": 0.1}, sharex = True)

for ibeam in ibeams:
    ax_top.scatter(df_all["runnum"], df_all[ibeam], s = 15, label=ibeam)
for ratio in ratios:
    x = df_all["runnum"]
    y = df_all[ratio].values

    p0 = p0_fit(y)

    scatter = ax_bot.scatter(x, y, s=15, label = ratio)
    color = scatter.get_facecolor()[0]

    ax_bot.axhline(p0, linestyle = "--", color = color, label = f"{ratio} avg: {100*p0:.2f}%")

ax_top.grid()
ax_top.legend()
ax_top.set_ylabel(f"Beam Currents $\mu A$")

ax_bot.grid()
ax_bot.legend()
ax_bot.set_ylabel(f"Beam Current Ratios")
ax_bot.set_xlabel("Run Number")

ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_bot.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

# plt.savefig(f"PNGs/bcms_vs_runnum.png", dpi=150, bbox_inches="tight")

# plt.show()

plt.close()

# ------------------------------------------------------
# Second plot, bcm current and ratios vs current
# ------------------------------------------------------
fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize = (8, 6), gridspec_kw = {"height_ratios":[1,4], "hspace": 0.1}, sharex = True)

for ibeam in ibeams:
    ax_top.scatter(df_all["ibeam2"], df_all[ibeam], s = 15, label=ibeam)
for ratio in ratios:
    x = df_all["ibeam2"]
    y = df_all[ratio].values

    p0 = p0_fit(y)

    scatter = ax_bot.scatter(x, y, s=15, label = ratio)
    color = scatter.get_facecolor()[0]

    ax_bot.axhline(p0, linestyle = "--", color = color, label = f"{ratio} avg: {100*p0:.2f}%")


ax_top.grid()
ax_top.legend()
ax_top.set_ylabel(f"Beam Currents $\mu A$")

ax_bot.grid()
ax_bot.legend()
ax_bot.set_ylabel(f"Beam Current Ratios")
ax_bot.set_xlabel("Beam Current (BCM2)")

ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_bot.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

# plt.savefig(f"PNGs/bcms_vs_current.png", dpi=150, bbox_inches="tight")

# plt.show()

plt.close()

# ------------------------------------------------------
# Third plot, looking at the different BCM2 calibrations
# ------------------------------------------------------
fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize = (8, 6), gridspec_kw = {"height_ratios":[1,4], "hspace": 0.1}, sharex = True)

for calib in calibs:
    ax_top.scatter(df_all["ibeam2"], df_all[calib], s = 15, label=calib)

ax_bot.scatter(df_all["ibeam2"], df_all["bcm2_calib_i1_ov_i2"], s = 15, label = "1 - calib_1_i / calib_2_i")

ax_top.grid()
ax_top.legend()
ax_top.set_ylabel(f"Beam Currents $\mu A$")

ax_bot.grid()
ax_bot.legend()
ax_bot.set_ylabel(f"1 - Beam Current Ratios")
ax_bot.set_xlabel("Beam Current (BCM2)")

ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax_bot.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
ax_bot.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))

# plt.savefig(f"PNGs/bcm2_calibs_vs_current.png", dpi=150, bbox_inches="tight")

plt.show()

plt.close()

# ------------------------------------------------------
# Forth plot, looking at the different BCM2 calibrations
# ------------------------------------------------------
fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize = (8, 6), gridspec_kw = {"height_ratios":[1,4], "hspace": 0.1}, sharex = True)

for calib in calibs:
    ax_top.scatter(df_all["runnum"], df_all[calib], s = 15, label=calib)

ax_bot.scatter(df_all["runnum"], df_all["bcm2_calib_i1_ov_i2"], s = 15, label = "1 - calib_1_i / calib_2_i")

ax_top.grid()
ax_top.legend()
ax_top.set_ylabel(f"Beam Currents $\mu A$")

ax_bot.grid()
ax_bot.legend()
ax_bot.set_ylabel(f"1 - Beam Current Ratios")
ax_bot.set_xlabel("Run Number")

ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax_bot.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
ax_bot.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))

# plt.savefig(f"PNGs/bcm2_calibs_vs_runnum.png", dpi=150, bbox_inches="tight")

plt.show()

plt.close()

# ------------------------------------------------------
# Looping over specific targets,
# ------------------------------------------------------

for target in targets:
    for beam_pass in beam_passes:
        df_targ = df_all[df_all["target"] == target]
        df_pass = df_targ[df_targ["beam_pass"] == beam_pass]
        # ------------------------------------------------------
        # Redoing the third plot now per target, looking at the different BCM2 calibrations
        # ------------------------------------------------------
        fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize = (8, 6), gridspec_kw = {"height_ratios":[1,4], "hspace": 0.1}, sharex = True)
    
        for calib in calibs:
            ax_top.scatter(df_pass["ibeam2"], df_pass[calib], s = 15, label=calib)

        ax_bot.scatter(df_pass["ibeam2"], df_pass["bcm2_calib_i1_ov_i2"], s = 15, label = "1 - calib_1_i / calib_2_i")

        ax_top.grid()
        ax_top.legend()
        ax_top.set_ylabel(f"Beam Currents $\mu A$")
        ax_top.set_title(f"{target.upper()} {beam_pass}: BCM2 Calibration Comparison")
        
        ax_bot.grid()
        ax_bot.legend()
        ax_bot.set_ylabel(f"1 - Beam Current Ratios")
        ax_bot.set_xlabel("Beam Current (BCM2)")
        
        ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        # ax_bot.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax_bot.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
        
        plt.savefig(f"PNGs/{target.lower()}_{beam_pass}_bcm2_calibs_vs_current.png", dpi=150, bbox_inches="tight")
    
        # plt.show()
    
        plt.close()

        # ------------------------------------------------------
        # Forth plot, looking at the different BCM2 calibrations
        # ------------------------------------------------------
        fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize = (8, 6), gridspec_kw = {"height_ratios":[1,4], "hspace": 0.1}, sharex = True)

        for calib in calibs:
            ax_top.scatter(df_pass["runnum"], df_pass[calib], s = 15, label=calib)

        ax_bot.scatter(df_pass["runnum"], df_pass["bcm2_calib_i1_ov_i2"], s = 15, label = "1 - calib_1_i / calib_2_i")

        ax_top.grid()
        ax_top.legend()
        ax_top.set_ylabel(f"Beam Currents $\mu A$")
        ax_top.set_title(f"{target.upper()} {beam_pass}: BCM2 Calibration Comparison")
    
        ax_bot.grid()
        ax_bot.legend()
        ax_bot.set_ylabel(f"1 - Beam Current Ratios")
        ax_bot.set_xlabel("Run Number")
    
        ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        # ax_bot.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax_bot.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))

        plt.savefig(f"PNGs/{target.lower()}_{beam_pass}_bcm2_calibs_vs_runnum.png", dpi=150, bbox_inches="tight")

        # plt.show()
    
        plt.close()
    



