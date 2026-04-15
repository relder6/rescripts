#!/usr/bin/env python3

import os
import csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
import mplhep
from collections import defaultdict

csv_files = ["CSVs/RME_Delta_R_results.csv", "CSVs/SLAC_E140_Delta_R.csv"]

all_data = []
for filepath in csv_files:
    if not os.path.exists(filepath):
        print(f"Skipping missing file: {filepath}")
        continue

    with open(filepath, "r", newline = "") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                all_data.append({"exp": row["exp"].strip(),
                                 "target": row["target"].strip().lower(),
                                 "xbj": float(row["xbj"]),
                                 "q2": float(row["q2_avg"]),
                                 "delta_R": float(row["delta_R"]),
                                 "delta_R_err": float(row["delta_R_err"]),
                                 "source": os.path.basename(filepath)})
            except (ValueError, KeyError, AttributeError):
                continue

xbj_groups = defaultdict(list)
for i, d in enumerate(all_data):
    xbj_groups[d["xbj"]].append(i)

xbj_offset = {}
step = 0.002

for xbj, idxs in xbj_groups.items():
    n = len(idxs)
    if n == 1:
        xbj_offset[idxs[0]] = 0.0
    else:
        start = -step * (n - 1) / 2.0
        for j, idx in enumerate(idxs):
            xbj_offset[idx] = start + j * step

exp_target_q2 = defaultdict(set)
for d in all_data:
    exp_target_q2[(d["exp"], d["target"])].add(round(d["q2"], 2))

has_multiple_q2 = {}
for (exp, tgt), q2s in exp_target_q2.items():
    has_multiple_q2[(exp, tgt)] = len(q2s) > 1

print("Loaded rows:", len(all_data))
if not all_data:
    raise SystemExit("No data found.")

plt.style.use(mplhep.style.ROOT)
plt.rcParams.update({"figure.titlesize": 26,
                     "axes.titlesize": 20,
                     "axes.labelsize": 18,
                     "legend.fontsize": 12,
                     "xtick.labelsize": 14,
                     "ytick.labelsize": 14})

fig, ax = plt.subplots(figsize = (5, 5))
# plt.subplots_adjust(right = 0.74, bottom = 0.18, left = 0.15)

base_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
target_list = sorted({d["target"] for d in all_data})
target_to_color = {tgt: base_colors[i % len(base_colors)] for i, tgt in enumerate(target_list)}
marker_cycle = ["o", "s", "^", "v", "D", "P", "X", "<", ">", "*", "h", "H", "p", "+", "x"]
target_to_marker = {tgt: marker_cycle[i % len(marker_cycle)] for i, tgt in enumerate(target_list)}

q2_min = min(d["q2"] for d in all_data)
q2_max = max(d["q2"] for d in all_data)

def q2_alpha(q2):
    if q2_max == q2_min:
        return 1.0
    return 0.5 + 0.5 * (q2 - q2_min) / (q2_max - q2_min)

def brighten(color, factor = 1.35):
    r, g, b, a = to_rgba(color)
    return (min(r * factor, 1), min(g * factor, 1), min(b * factor, 1), a)

seen_labels = set()
legend_handles = []

for i, d in enumerate(all_data):
    tgt = d["target"]
    exp = d["exp"]
    q2 = d["q2"]
    is_rme = (d["source"] == "RME_Delta_R_results.csv")

    base_color = target_to_color[tgt]
    key = (d["exp"], d["target"])
    if has_multiple_q2[key]:
        alpha = q2_alpha(q2)
    else:
        alpha = 1.0
    edge_color = brighten(base_color) if is_rme else base_color
    face_color = brighten(base_color) if is_rme else "none"
    size = 130 if is_rme else 70
    zorder = 3 if is_rme else 2

    label = f"{exp} {tgt} Q2 = {q2:.2f}"

    x_plot = d["xbj"] + xbj_offset[i]

    ax.errorbar(x_plot, d["delta_R"],yerr = d["delta_R_err"],fmt = "none",ecolor = (*to_rgba(base_color)[:3], alpha), capsize = 0,zorder = 1,)

    ax.scatter(x_plot, d["delta_R"],marker = target_to_marker[tgt],s = size,
               facecolors = face_color if face_color == "none" else (*to_rgba(face_color)[:3], alpha),edgecolors = ("black" if is_rme else (*to_rgba(edge_color)[:3], alpha)),linewidths = 1.5,zorder = zorder,)

    
    
    if label not in seen_labels:
        seen_labels.add(label)
        legend_handles.append(Line2D([0], [0],marker = target_to_marker[tgt],linestyle = "none",markerfacecolor = face_color if face_color == "none" else (*to_rgba(face_color)[:3], alpha),markeredgecolor = ("black" if is_rme else (*to_rgba(edge_color)[:3], alpha)),color = "w",label = label,markersize = 8 if is_rme else 6))

ax.axhline(y = 0, color = "black", linewidth = 2, zorder = 0)
ax.set_xlabel(r"$x_{bj}$")
ax.set_ylabel(r"$\Delta R = R_A - R_D$")
ax.set_title(r"Comparison of $\Delta R = R_A - R_D$")
ax.grid()

xmin = min(p["xbj"] for p in all_data)
xmax = max(p["xbj"] for p in all_data)
ymin = min(p["delta_R"] - p["delta_R_err"] for p in all_data)
ymax = max(p["delta_R"] + p["delta_R_err"] for p in all_data)
pad_frac = 0.05
# # ax.set_xlim(xmin - pad_frac * (xmax - xmin), xmax + pad_frac * (xmax - xmin))
# ax.set_xlim(0,  xmax + pad_frac * (xmax - xmin))
# ax.set_ylim(ymin - pad_frac * (ymax - ymin), ymax + pad_frac * (ymax - ymin))

ax.legend(handles = legend_handles, title = "Legend", loc = "best", title_fontsize = 14, frameon = True, fancybox = True,
          framealpha = 0.9, facecolor = "white", edgecolor = "black")
plt.tight_layout()
plt.show()
