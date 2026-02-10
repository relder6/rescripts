#!/usr/bin/env python3

import os, sys
import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from datetime import datetime

# ----------------------------------------------
# CSV files
# ----------------------------------------------
csv_files = ["CSVs/RME_emc_results.csv", "CSVs/WorldData.csv", "CSVs/DG_emc.csv"]

# ----------------------------------------------
# Target definitions, mapping, selection
# ----------------------------------------------
targets = {
    "c":   {"A": 12,  "Z": 6,  "name": "Carbon",     "symbol": "C",   "aliases": ["c12","carbon"]},
    "al":  {"A": 27,  "Z": 13, "name": "Aluminum",   "symbol": "Al",  "aliases": ["al13","aluminum"]},
    "cu":  {"A": 64,  "Z": 29, "name": "Copper",     "symbol": "Cu",  "aliases": ["cu29","copper"]},
    "fe":  {"A": 56,  "Z": 26, "name": "Iron",       "symbol": "Fe",  "aliases": ["iron"]},
    "pb":  {"A": 207, "Z": 82, "name": "Lead",       "symbol": "Pb",  "aliases": ["lead"]},
    "ld2": {"A": 2,   "Z": 1,  "name": "Deuterium",  "symbol": "D",   "aliases": ["d2","ld2"]},
    "lh2": {"A": 1,   "Z": 1,  "name": "Hydrogen",   "symbol": "H",   "aliases": ["h2","lh2"]},
}

alias_to_key = {}
for key, info in targets.items():
    alias_to_key[key] = key
    for a in info.get("aliases", []):
        alias_to_key[a] = key

if len(sys.argv) == 2:
    selected_input = sys.argv[1].strip().lower()
else:
    selected_input = input("Enter desired target for comparison: ").strip().lower()

selected_key = alias_to_key.get(selected_input)
if not selected_key:
    print(f"Unknown target: {selected_input}")
    exit(1)

target_info = targets[selected_key]
A = target_info["A"]
Z = target_info["Z"]
A_nominal = int(round(A))
Z_int = int(round(Z))
target_name = target_info["name"]
target_symbol = target_info["symbol"]

# ----------------------------------------------
# Handling CSV inputs
# ----------------------------------------------
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
                    "exp": row["exp"].strip().upper(),
                    "A": float(row["A"]),
                    "Z": float(row["Z"]),
                    "xbj": float(row["xbj"]),
                    "emc_ratio": float(row["emc_ratio"]),
                    "emc_ratio_err": float(row["emc_ratio_err"]),
                    "source": os.path.basename(filepath)
                })
            except (ValueError, KeyError):
                continue

filtered_data = [d for d in all_data if int(round(d["A"])) == A_nominal and int(round(d["Z"])) == Z_int]

if not filtered_data:
                 print(f"No data found for A={A_nominal}, Z={Z_int}")
                 exit(1)
                 
experiments = {}
for d in filtered_data:
    experiments.setdefault(d["exp"], []).append(d)

unique_sources = sorted({d["source"] for d in filtered_data})

# ----------------------------------------------
# Plotting
# ----------------------------------------------
fig, ax = plt.subplots(figsize=(7,5))
plt.subplots_adjust(bottom=0.18, left=0.18)

marker_cycle = ['o', 's', '^', 'v', 'P', 'X', '*', '<', '>']
source_to_marker = {src: ("D" if src=="RME_results.csv" else marker_cycle.pop(0)) for src in unique_sources}

color_cycle = ["#000000","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf"]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_cycle)
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

lines = {}
for exp, points_all in experiments.items():
    by_source = {}
    for p in points_all:
        by_source.setdefault(p["source"], []).append(p)

    for src, points in by_source.items():
        color = color_cycle[len(lines) % len(color_cycle)]
        lines[(exp, src)] = True
        x = np.array([p["xbj"] for p in points])
        y = np.array([p["emc_ratio"] for p in points])
        yerr = np.array([p["emc_ratio_err"] for p in points])
        order = np.argsort(x)
        x, y, yerr = x[order], y[order], yerr[order]
        ms = 4 if src=="RME_emc_results.csv" else 3
        ax.errorbar(x, y, yerr=yerr, fmt=source_to_marker[src], linestyle='none',
                    markersize=ms, capsize=0, markerfacecolor='none', markeredgewidth=1.1,
                    color=color, label=exp)

ax.set_xlabel(r"$x_{bj}$", fontsize = 14)
ax.set_ylabel(rf"$(\sigma_{{{target_symbol}}} \, / \, {A_nominal}) \: / \: (\sigma_D \, / \, 2)$", fontsize=14)
ax.legend()
ax.grid(alpha=0.3)

xmin, xmax = min(p["xbj"] for p in filtered_data), max(p["xbj"] for p in filtered_data)
ymin, ymax = min(p["emc_ratio"] - p["emc_ratio_err"] for p in filtered_data), max(p["emc_ratio"] + p["emc_ratio_err"] for p in filtered_data)
pad_frac = 0.05
ax.set_xlim(xmin - pad_frac*(xmax-xmin), xmax + pad_frac*(xmax-xmin))
ax.set_ylim(ymin - pad_frac*(ymax-ymin), ymax + pad_frac*(ymax-ymin))

# ----------------------------------------------
# Sliders
# ----------------------------------------------
h, gap = 0.015, 0.01
ax_xmin = plt.axes([0.25, 0.08 + h + gap, 0.65, h])
ax_xmax = plt.axes([0.25, 0.08, 0.65, h])
ax_ymin = plt.axes([0.08, 0.25, h, 0.65])
ax_ymax = plt.axes([0.11, 0.25, h, 0.65])

slider_xmin = Slider(ax_xmin, "Xmin", xmin, xmax, valinit=xmin)
slider_xmax = Slider(ax_xmax, "Xmax", xmin, xmax, valinit=xmax)
slider_ymin = Slider(ax_ymin, "Ymin", ymin, ymax, valinit=ymin, orientation="vertical")
slider_ymax = Slider(ax_ymax, "Ymax", ymin, ymax, valinit=ymax, orientation="vertical")

for s in [slider_xmin, slider_xmax, slider_ymin, slider_ymax]:
    s.poly.set_alpha(0.25)
    s.valtext.set_visible(False)
    s.ax.set_facecolor("0.95")

def update(val):
    xmin, xmax = slider_xmin.val, slider_xmax.val
    ymin, ymax = slider_ymin.val, slider_ymax.val
    if xmin >= xmax or ymin >= ymax: return
    ax.set_xlim(xmin - pad_frac*(xmax-xmin), xmax + pad_frac*(xmax-xmin))
    ax.set_ylim(ymin - pad_frac*(ymax-ymin), ymax + pad_frac*(ymax-ymin))
    fig.canvas.draw_idle()

for s in [slider_xmin, slider_xmax, slider_ymin, slider_ymax]:
    s.on_changed(update)

# ----------------------------------------------
# Save button
# ----------------------------------------------
ax_save = plt.axes([0.01, 0.01, 0.1, 0.05])
btn_save = Button(ax_save, "Save PNG", color='lightgray', hovercolor='0.975')

def save_png(event):
    for s in [slider_xmin, slider_xmax, slider_ymin, slider_ymax]:
        s.ax.set_visible(False)
    ax_save.set_visible(False)
    os.makedirs("PNGs", exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    save_path = os.path.join("PNGs", f"EMC_{selected_key}_{timestamp}.png")
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    for s in [slider_xmin, slider_xmax, slider_ymin, slider_ymax]:
        s.ax.set_visible(True)
    ax_save.set_visible(True)
    print(f"Saved figure as {save_path}")

btn_save.on_clicked(save_png)

plt.show()
