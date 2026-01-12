#!/usr/bin/env python3
import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# ----------------------
# CSV files to load
# ----------------------
csv_files = ["WorldData.csv", "RME_results.csv", "DG_carbon.csv"]

# ----------------------
# Targets dictionary
# ----------------------
targets = {
    "carbon": (12, 6), "c": (12, 6),
    "aluminum": (27, 13), "al": (27, 13),
    "iron": (56, 26), "fe": (56, 26),
    "lead": (208, 82), "pb": (208, 82),
    "copper": (64, 29), "cu": (64, 29)
}

# ----------------------
# Choose target
# ----------------------
if len(sys.argv) == 2:
    selected_target = sys.argv[1].strip().lower()
else:
    selected_target = input("Enter desired target for comparison: ").strip().lower()

if selected_target not in targets:
    print(f"Target {selected_target} unknown or not supported.")
    exit(1)

A, Z = targets[selected_target]

selected_target_shortname_to_title_longname = {
    "al": "Aluminum",
    "c": "Carbon",
    "cu": "Copper",
    "opt1": "Optics1",
    "opt2": "Optics2",
    "ld2": "Deuterium",
    "lh2": "Hydrogen",
    "hole": "Carbon Hole",
    "dummy": "Dummy"
}
selected_target_titlename = selected_target_shortname_to_title_longname.get(selected_target, selected_target)

# ----------------------
# Load data
# ----------------------
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
                    "source": os.path.basename(filepath)
                })
            except ValueError:
                continue

# ----------------------
# Filter data for selected target
# ----------------------
filtered_data = [d for d in all_data if abs(d["A"]-A)<1e-6 and abs(d["Z"]-Z)<1e-6]
if len(filtered_data) == 0:
    print(f"No world data found for A={A}, Z={Z}")
    exit(1)

# ----------------------
# Setup experiments and unique CSV sources
# ----------------------
experiments = {}
for data in filtered_data:
    experiments.setdefault(data["exp"], []).append(data)

marker_cycle = ['o', 's', '^', 'v', 'D', 'P', 'X', '*', '<', '>']
unique_sources = sorted({d["source"] for d in filtered_data})
source_to_marker = {src: marker_cycle[i % len(marker_cycle)] for i, src in enumerate(unique_sources)}

color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
exp_list = sorted(experiments.keys())
exp_to_color = {exp: color_cycle[i % len(color_cycle)] for i, exp in enumerate(exp_list)}

# ----------------------
# Figure and axes
# ----------------------
fig, ax = plt.subplots(figsize=(7, 5))
plt.subplots_adjust(bottom=0.3, left=0.18)

ax.set_xlabel(r"$x_{Bj}$")
ax.set_ylabel("Ratio")
ax.set_title(f"{selected_target_titlename} World Data Comparison")
ax.grid(alpha=0.3)

xmin_init = np.min([p["xbj"] for p in filtered_data])
xmax_init = np.max([p["xbj"] for p in filtered_data])
ymin_init = np.min([p["ratio"] - p["TotErrAbs"] for p in filtered_data])
ymax_init = np.max([p["ratio"] + p["TotErrAbs"] for p in filtered_data])
padding_frac = 0.05
x_pad = padding_frac * (xmax_init - xmin_init)
y_pad = padding_frac * (ymax_init - ymin_init)
ax.set_xlim(xmin_init - x_pad, xmax_init + x_pad)
ax.set_ylim(ymin_init - y_pad, ymax_init + y_pad)

# ----------------------
# Sliders
# ----------------------
h = 0.015
gap = 0.01

ax_xmin = plt.axes([0.25, 0.08 + 3*(h+gap), 0.65, h])
ax_xmax = plt.axes([0.25, 0.08 + 2*(h+gap), 0.65, h])
ax_scale = plt.axes([0.25, 0.08 + 1*(h+gap), 0.65, h])
ax_ymin = plt.axes([0.08, 0.25, h, 0.65])
ax_ymax = plt.axes([0.11, 0.25, h, 0.65])

slider_xmin = Slider(ax_xmin, "Xmin", xmin_init, xmax_init, valinit=xmin_init)
slider_xmax = Slider(ax_xmax, "Xmax", xmin_init, xmax_init, valinit=xmax_init)
slider_scale = Slider(ax_scale, "RME Scale", -0.1, 0.1, valinit=0.0)
slider_ymin = Slider(ax_ymin, "Ymin", ymin_init, ymax_init, valinit=ymin_init, orientation="vertical")
slider_ymax = Slider(ax_ymax, "Ymax", ymin_init, ymax_init, valinit=ymax_init, orientation="vertical")

for s in [slider_xmin, slider_xmax, slider_ymin, slider_ymax]:
    s.poly.set_alpha(0.25)
    s.valtext.set_visible(False)
    s.ax.set_facecolor("0.95")
  
slider_scale.poly.set_alpha(0.25)
slider_scale.ax.set_facecolor("0.95")

# ----------------------
# Update function
# ----------------------
def update(val):
    xmin, xmax = slider_xmin.val, slider_xmax.val
    ymin, ymax = slider_ymin.val, slider_ymax.val
    scale = slider_scale.val

    slider_scale.valtext.set_text(f"{scale*100:+0.1f}%")

    if xmin >= xmax or ymin >= ymax:
        return

    xpad = padding_frac * (xmax - xmin)
    ypad = padding_frac * (ymax - ymin)
    ax.set_xlim(xmin - xpad, xmax + xpad)
    ax.set_ylim(ymin - ypad, ymax + ypad)

    # Clear previous lines
    ax.collections.clear()
    ax.lines.clear()

    if ax.legend_:
        ax.legend_.remove()

    added_labels = set()

    # Draw all experiments
    for exp in exp_list:
        points_all = experiments[exp]
        by_source = {}
        for p in points_all:
            by_source.setdefault(p["source"], []).append(p)

        for src, points in by_source.items():
            x = np.array([p["xbj"] for p in points])
            y = np.array([p["ratio"] for p in points])
            yerr = np.array([p["TotErrAbs"] for p in points])

            if src == "RME_results.csv":
                y = y * (1 + scale)
                yerr = yerr * (1 + scale)

            order = np.argsort(x)
            x, y, yerr = x[order], y[order], yerr[order]

            ax.errorbar(
                x, y, yerr=yerr,
                fmt=source_to_marker[src],
                color=exp_to_color[exp],
                linestyle='none',
                markersize=4,
                capsize=0,
                markerfacecolor='none',
                markeredgewidth=1.1,
                label=exp if src == list(by_source.keys())[0] else None
            )

    fig.canvas.draw_idle()

# ----------------------
# Connect sliders
# ----------------------
slider_xmin.on_changed(update)
slider_xmax.on_changed(update)
slider_ymin.on_changed(update)
slider_ymax.on_changed(update)
slider_scale.on_changed(update)

# ----------------------
# Initial draw
# ----------------------
update(None)
plt.show()
