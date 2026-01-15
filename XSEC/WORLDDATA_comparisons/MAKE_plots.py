#!/usr/bin/env python3

import os, re, sys
import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# csv_files = ["WorldData.csv", "RME_results.csv", "DG_carbon.csv", "SCALED_RME_results.csv"]
# csv_files = ["WorldData.csv", "SCALED_RME_results.csv"]
csv_files = ["WorldData.csv", "RME_results.csv"]

targets = {
    "carbon": (12, 6),
    "c": (12, 6),
    "aluminum": (27, 13),
    "al": (27, 13),
    "iron": (56, 26),
    "fe": (56, 26),
    "lead": (208, 82),
    "pb": (208, 82),
    "copper": (64, 29),
    "cu": (64, 29)
}

if len(sys.argv) == 2:
    selected_target = sys.argv[1].strip().lower()
else:
    selected_target = input("Enter desired target for comparison: ").strip().lower()

if selected_target not in targets:
    print(f"Target {selected_target} unknown or not currently supported; please try again.")
    exit(1)

selected_target_shortcut_to_target_variable = {"al":"al","al13":"al","aluminum":"al",
                                                   "c":"c","c12":"c","carbon":"c",
                                                   "cu":"cu","cu29":"cu","copper":"cu",
                                                   "opt1":"optics1","optics1":"optics1",
                                                   "opt2":"optics2","optics2":"optics2",
                                                   "d2":"ld2","ld2":"ld2",
                                                   "h2":"lh2","lh2":"lh2",
                                                   "hole":"hole","chole":"hole","c-hole":"hole",
                                                   "dummy":"dummy","dum":"dummy",
                                                   }

selected_target_shortname = selected_target_shortcut_to_target_variable.get(selected_target)

if not selected_target_shortname:
        print(f"Unknown target: {selected_target}.  Please try again.")
        exit(1)

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

selected_target_titlename = selected_target_shortname_to_title_longname.get(selected_target_shortname)

A, Z = targets[selected_target]

all_data = []

for filepath in csv_files:
    if not os.path.exists(filepath):
        print(f"Skipping missing file: {filepath}")
        continue
    
    with open(filepath, "r") as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            try:
                all_data.append({
                    "exp": row["Exp"].strip().upper(),
                    "A": float(row["A"]),
                    "Z": float(row["Z"]),
                    "xbj": float(row["x"]),
                    "ratio": float(row["Ratio"]),
                    "TotErrAbs": float(row["TotalErrorAbs"]),
                    "source": os.path.basename(filepath)})
                    # "xi": float(row["xi"]),
                    # "Ebeam": float(row["beamE(GeV)"]),
                    # "theta": float(row["angle"]),
                    # "StatAbs": float(row["Stat(abs)"]),
                    # "SysAbs": float(row["Sys(abs)"]),
                    # "Sys_xcorr": float(row["Sys-xcorr"]),
                    # "iso": float(row["Iso"]),
                    # "coul": float(row["Coul"]),
                    # "ratio_unc": float(row["UncRatio"]),
                    # "NormUncPct": float(row["NormUncPct"]),
            except ValueError:
                continue

filtered_data = []

for data in all_data:
    if abs(data["A"] - A) < 1e-6 and abs(data["Z"] - Z) < 1e-6:
        filtered_data.append(data)

if len(filtered_data) == 0:
    print(f"No world data found for A={A}, Z={Z}")
    exit(1)

experiments = {}

for data in filtered_data:
    exp = data["exp"]
    if exp not in experiments:
        experiments[exp] = []

    experiments[exp].append(data)
    
unique_sources = sorted({d["source"] for d in filtered_data})



# ----------------------
# Making the figure here
# ----------------------
fig, ax = plt.subplots(figsize=(7, 5))
plt.subplots_adjust(bottom=0.18, left=0.18)

marker_cycle = ['o', 's', '^', 'v', 'D', 'P', 'X', '*', '<', '>']

source_to_marker = {src: marker_cycle[i % len(marker_cycle)] for i, src in enumerate(unique_sources)}

plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.tab10.colors)

lines = {}

for exp, points_all in experiments.items():

    by_source = {}
    for p in points_all:
        by_source.setdefault(p["source"], []).append(p)

    for src, points in by_source.items():
        x = np.array([p["xbj"] for p in points])
        y = np.array([p["ratio"] for p in points])
        yerr = np.array([p["TotErrAbs"] for p in points])

        order = np.argsort(x)
        x, y, yerr = x[order], y[order], yerr[order]

        ax.errorbar(x, y, yerr=yerr, fmt=source_to_marker[src], linestyle='none', markersize=4, capsize=0, markerfacecolor='none', markeredgewidth=1.1, label=exp)

ax.set_xlabel(r"$x_{Bj}$")
ax.set_ylabel("Ratio")
ax.set_title(f"{selected_target_titlename} World Data Comparison")
ax.legend()
ax.grid(alpha=0.3)

xmin_init = np.min([p["xbj"] for p in filtered_data])
xmax_init = np.max([p["xbj"] for p in filtered_data])
ymin_init = np.min([p["ratio"] - p["TotErrAbs"] for p in filtered_data])
ymax_init = np.max([p["ratio"] + p["TotErrAbs"] for p in filtered_data])

padding_frac = 0.05

x_pad = padding_frac * (xmax_init - xmin_init)
y_pad = padding_frac * (ymax_init - ymin_init)

xmin_pad = xmin_init - x_pad
xmax_pad = xmax_init + x_pad
ymin_pad = ymin_init - y_pad
ymax_pad = ymax_init + y_pad

ax.set_xlim(xmin_init - x_pad, xmax_init + x_pad)
ax.set_ylim(ymin_init - y_pad, ymax_init + y_pad)

# Adding the slider bars
h = 0.015
gap = 0.01

ax_xmin = plt.axes([0.25, 0.08 + h + gap, 0.65, h])
ax_xmax = plt.axes([0.25, 0.08,           0.65, h])

ax_ymin = plt.axes([0.08, 0.25, h, 0.65])
ax_ymax = plt.axes([0.11, 0.25, h, 0.65])

slider_xmin = Slider(ax_xmin, "Xmin", xmin_init, xmax_init, valinit=xmin_init)
slider_xmax = Slider(ax_xmax, "Xmax", xmin_init, xmax_init, valinit=xmax_init)
slider_ymin = Slider(ax_ymin, "Ymin", ymin_init, ymax_init, valinit=ymin_init, orientation="vertical")
slider_ymax = Slider(ax_ymax, "Ymax", ymin_init, ymax_init, valinit=ymax_init, orientation="vertical")

for s in [slider_xmin, slider_xmax, slider_ymin, slider_ymax]:
    s.poly.set_alpha(0.25)
    s.valtext.set_visible(False)
    s.ax.set_facecolor("0.95")
    
def update(val):
    xmin = slider_xmin.val
    xmax = slider_xmax.val
    ymin = slider_ymin.val
    ymax = slider_ymax.val

    if xmin >= xmax or ymin >= ymax:
        return

    xpad = padding_frac * (xmax - xmin)
    ypad = padding_frac * (ymax - ymin)

    ax.set_xlim(xmin - xpad, xmax + xpad)
    ax.set_ylim(ymin - ypad, ymax + ypad)
    fig.canvas.draw_idle()

slider_xmin.on_changed(update)
slider_xmax.on_changed(update)
slider_ymin.on_changed(update)
slider_ymax.on_changed(update)

plt.show()
                
