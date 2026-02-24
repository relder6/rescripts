#!/usr/bin/env python3

import sys, os
import matplotlib
matplotlib.use('Agg')
import uproot
import numpy as np
import matplotlib.pyplot as plt
import boost_histogram as bh
import matplotlib.colors as mplcolors
from scipy.optimize import curve_fit

# --------------------------------------------------------------------------
# Check file locations before running,
# --------------------------------------------------------------------------
root_directory = f"/work/hallc/c-rsidis/relder/hallc_replay_rsidis/ROOTfiles"

# --------------------------------------------------------------------------
# Handling user inputs, unloading branches into numpy arrays
# --------------------------------------------------------------------------
if len(sys.argv) > 1:
    runnum = sys.argv[1]
else:
    runnum = input("Input the run number you wish to analyze: ")

coin_pattern = f"coin_replay_production_{runnum}_50000.root"
shms_pattern = f"shms_coin_replay_production_{runnum}_-1.root"
hms_pattern = f"hms_coin_replay_production_{runnum}_-1.root"

coin_path = f"{root_directory}/{coin_pattern}"
shms_path = f"{root_directory}/{shms_pattern}"
hms_path = f"{root_directory}/{hms_pattern}"

file_type = None
if os.path.exists(coin_path):
    input_root_filepath = coin_path
    file_type = "COIN"
    os.makedirs(file_type, exist_ok=True)
elif os.path.exists(shms_path):
    input_root_filepath = shms_path
    file_type = "SHMS"
    os.makedirs(file_type, exist_ok=True)
elif os.path.exists(hms_path):
    print(f"WARNING:\tRun {runnum} is an HMS replay; skipping...")
    sys.exit(1)
else:
    print(f"ERROR:\tRun {runnum} is missing a ROOTfile; skipping...")
    sys.exit(1)

d_calo_fp = 292.64 # distance from focal plane to calorimeter face

branches = ["P.dc.x_fp", "P.dc.y_fp", "P.dc.xp_fp", "P.dc.yp_fp", "P.gtr.dp", "P.hgcer.npeSum", "P.ngcer.npeSum", "P.cal.etottracknorm", "P.gtr.beta"]

data_cut = ("(P.gtr.dp > -15) & (P.gtr.dp < 27) & (P.hgcer.npeSum > 1.5) & (P.ngcer.npeSum > 1.5) & (P.gtr.beta > 0.8) & (P.gtr.beta < 1.2) & (P.cal.etottracknorm > 0)")

try:
    with uproot.open(input_root_filepath) as rootfile:
        tree = rootfile["T"]
        arrays = tree.arrays(branches, cut=data_cut, library="np")
except uproot.exceptions.KeyInFileError as e:
    print(f"ERROR:\tRun {runnum} is a {file_type} file missing branch '{e.key}'; skipping...")
    sys.exit(1) 
    
# --------------------------------------------------------------------------
# Setting up variables for plotting
# --------------------------------------------------------------------------

xcalo = (arrays["P.dc.x_fp"]) + (arrays["P.dc.xp_fp"])*d_calo_fp

ycalo = (arrays["P.dc.y_fp"]) + (arrays["P.dc.yp_fp"])*d_calo_fp

weight = arrays["P.cal.etottracknorm"]

finite_mask = np.isfinite(xcalo) & np.isfinite(ycalo) & np.isfinite(weight)

xcalo = xcalo[finite_mask]

ycalo = ycalo[finite_mask]

weight = weight[finite_mask]

xbins, ybins = 100, 100

xmin, xmax = -60, 60

ymin, ymax = -58, 58

xrange = (xmin - 10, xmax + 10)

yrange = (ymin - 10, ymax + 10)

hCaloPos = bh.Histogram(
    bh.axis.Regular(xbins, xrange[0], xrange[1]),
    bh.axis.Regular(ybins, yrange[0], yrange[1]))

hCaloPosWt = bh.Histogram(
    bh.axis.Regular(xbins, xrange[0], xrange[1]),
    bh.axis.Regular(ybins, yrange[0], yrange[1]),
    storage = bh.storage.Weight())

hCaloPos.fill(ycalo, xcalo)
hCaloPosWt.fill(ycalo, xcalo, weight=weight)

counts_unweighted = hCaloPos.values()
counts_weighted = hCaloPosWt.values()

total_counts = counts_unweighted.sum()
if total_counts == 0:
    print(f"WARNING:\tRun {runnum} has zero tracks; skipping...")
    sys.exit(2)

with np.errstate(divide='ignore', invalid='ignore'):
    hCaloPosNormU = np.divide(counts_weighted, counts_unweighted)
    hCaloPosNormU[~np.isfinite(hCaloPosNormU)] = 0  # set inf/NaN to 0
    
# -----------------------------------------------------------------------------
# Plotting figures, disabling them  but will leave in code for debugging
# -----------------------------------------------------------------------------
# plt.figure(figsize=(8, 6))
# plt.imshow(
#     hCaloPos.view().T,
#     origin="lower",
#     extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
#     aspect="auto"
# )
# plt.xlabel("xcalo")
# plt.ylabel("ycalo")
# plt.title("2D Histogram: ycalo vs xcalo")
# plt.colorbar(label="Counts")
# plt.savefig("hCaloPos.png", bbox_inches="tight", dpi=150)
# plt.close()

# plt.figure(figsize=(8, 6))
# plt.imshow(
#     hCaloPosWtU.view().T,
#     origin="lower",
#     extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
#     aspect="auto"
# )
# plt.xlabel("xweight")
# plt.ylabel("yweight")
# plt.title("yweight v xweight")
# plt.colorbar(label="Counts")
# plt.savefig("hCaloPosWtU.png", bbox_inches="tight", dpi=150)
# plt.close()

# -----------------------------------------------------------------------------
# Getting the figure ready to place both plots
# -----------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (14, 6), constrained_layout = True)
fig.suptitle(f"{file_type} Run {runnum}", fontsize=16, fontweight="bold")

# -----------------------------------------------------------------------------
# Plotting normalized e/p per track at calorimeter
# -----------------------------------------------------------------------------
data = hCaloPosNormU
vmin = 0.0
vmax_candidate = np.nanmax(data)
vmax = vmax_candidate if vmax_candidate > vmin else vmin + 1  # avoiding zero range here

def safe_norm_pos(val, vmin, vmax, prev_pos):
    if vmax == vmin:
        return 1.0
    pos = (val - vmin) / (vmax - vmin)
    pos = np.clip(pos, 0.0, 1.0)
    return max(pos, prev_pos)

pos_0 = 0.0
pos_1 = safe_norm_pos(1, vmin, vmax, pos_0)
pos_2 = safe_norm_pos(2, vmin, vmax, pos_1)
pos_max = 1.0

cmap = mplcolors.LinearSegmentedColormap.from_list(
    "custom_cmap",
    [
        (pos_0, "white"),
        (pos_1, "navy"),
        (pos_2, "orange"),
        (pos_max, "red")
    ]
)

norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)

im = ax1.imshow(
    data.T,
    origin="lower",
    extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
    aspect="auto",
    cmap=cmap,
    norm=norm
)

cbar = fig.colorbar(im, ax=ax1, format='%.1f', pad=0.01)

ticks = np.arange(0, int(np.floor(vmax)) + 0.5, 1)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])
cbar.ax.tick_params(labelsize=8)

ax1.set_title(f"Normalized E/p per Track at Calorimeter", fontsize = 14)

# Now adding the block grid,
numrows, numcols, blockspacing = 16, 14, 9
startrow, startcol = -blockspacing * numrows / 2.0, -blockspacing * numcols / 2.0
for i in range(numrows + 1):
    y = startrow + i * blockspacing
    ax1.hlines(y, startcol, startcol + blockspacing * numcols, colors='silver', linewidth=1.0, alpha=0.6)
for i in range(numcols + 1):
    x = startcol + i * blockspacing
    ax1.vlines(x, startrow, startrow + blockspacing * numrows, colors='silver', linewidth=1.0, alpha=0.6)

# -----------------------------------------------------------------------------
# Plotting fitted distribution of e/p
# -----------------------------------------------------------------------------
bin_min, bin_max, bin_num = 0, 2, 100
data_bins = np.linspace(bin_min, bin_max, bin_num + 1)
counts, bin_edges = np.histogram(arrays["P.cal.etottracknorm"], bins=data_bins)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Initializing fit_failed flag, several checks will be made against this to avoid wasting time fitting nothing
fit_failed = False

# Need a broad search window here; want to capture peaks that are poorly fit (--> bad calibration!)
search_min, search_max = 0.80, 1.20
search_mask = (bin_centers >= search_min) & (bin_centers <= search_max)

x_search = bin_centers[search_mask]
y_search = counts[search_mask]

nonzero = y_search > 0

x_search = x_search[nonzero]
y_search = y_search[nonzero]

if x_search.size == 0:
    print(f"WARNING:\tRun {runnum} no entries in search window; skipping...")
    fit_failed = True

# Here I'm finding the peak in the broad search window, and only keeping 0.15 around that,

if not fit_failed:    
    peak_idx = np.argmax(y_search)
    peak_x = x_search[peak_idx]
    peak_y = y_search[peak_idx]

    fit_width = 0.15
    local_mask = np.abs(x_search - peak_x) <= fit_width

    x_fit = x_search[local_mask]
    y_fit = y_search[local_mask]

    if x_fit.size > 0:
        fit_bin_min = np.min(y_fit)
        fit_bin_max = np.max(y_fit)
        fit_bin_avg = np.mean(y_fit)
        fit_bin_sum = np.sum(y_fit)
    else:
        fit_bin_min = fit_bin_max = fit_bin_avg = fit_bin_sum = np.nan

    if y_fit.size == 0:
        print(f"WARNING:\tRun {runnum} contains no events in range; skipping...")
        fit_failed = True

# This helps get rid of fitting random noise,
min_peak_height = 15
if peak_y < min_peak_height:
    print(f"WARNING: Run {runnum} peak too small ({peak_y}); skipping fit...")
    fit_failed = True

# We want to fit the Gaussian only around that found peak, keep_above cuts values /less than/ keep_above * peak_height
peak_height = np.max(y_fit)
keep_above = 0.30
top_mask = y_fit >= keep_above * peak_height

x_fit= x_fit[top_mask]
y_fit= y_fit[top_mask]

min_points = 3
if x_fit.size < min_points:
    print(f"WARNING:\tRun {runnum} not enough bins around peak top; skipping fit...")
    fit_failed = True
    
amp_guess = np.max(y_fit)
mean_guess = np.sum(x_fit * y_fit) / np.sum(y_fit)
sigma_guess = np.sqrt(np.sum(y_fit * (x_fit - mean_guess) ** 2) / np.sum(y_fit))

# min_peak_bin = 10

# if np.max(y_fit) < min_peak_bin:
#     print(f"WARNING:\tRun {runnum} peak bin too low ({np.max(y_fit)}); skipping fit...")
#     fit_failed = True

# low_peak_bin = 40

# if min_peak_bin < np.max(y_fit) < low_peak_bin:
#     bin_ratio = np.max(y_fit) / np.min(y_fit)
#     if bin_ratio < 1.3:
#         print(f"WARNING:\tRun {runnum} low peak, does not seem significant within region; skipping fit...")
#         fit_failed = True

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

if not fit_failed:
    try:
        popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=[amp_guess, mean_guess, sigma_guess])
        amp_fit, mean_fit, sigma_fit = popt
        amp_err, mean_err, sigma_err = np.sqrt(np.diag(pcov))
    except Exception:
        fit_failed = True
        print(f"WARNING:\tRun {runnum} was unable to find a Gaussian fit.")
        amp_fit = mean_fit = sigma_fit = amp_err = mean_err = sigma_err = np.nan
else:
    amp_fit = mean_fit = sigma_fit = amp_err = mean_err = sigma_err = np.nan

if not fit_failed:
    if not (search_min <= mean_fit <= search_max):
        fit_failed = True
        print(f"WARNING:\tRun {runnum} fitted mean {mean_fit:.3f} outside search window; skipping fit.")
        amp_fit = mean_fit = sigma_fit = amp_err = mean_err = sigma_err = np.nan
        
# Now plotting to the second pad
ax2.hist(arrays["P.cal.etottracknorm"], bins=data_bins, histtype='step', color='red', label='E/p')

if not fit_failed:
    x_plot = np.linspace(np.min(x_fit), np.max(x_fit), 500)
    ax2.plot(x_plot, gaussian(x_plot, *popt), color='navy', label='Gaussian fit')
    ax2.text(0.05, 0.95,
             rf"$\mu_{{\rm fit}} = {mean_fit:.3f} \pm {mean_err:.3f}$" + "\n" +
             rf"$\sigma_{{\rm fit}} = {sigma_fit:.3f} \pm {sigma_err:.3f}$",
             # rf"$A_{{\rm fit}} = {amp_fit:.3g} \pm {amp_err:.3g}$",
             transform=ax2.transAxes, verticalalignment='top', fontsize=10,
             bbox=dict(facecolor='white', alpha=0.7, edgecolor='black'))
else:
    ax2.text(0.5, 0.5, "FIT FAILED", color="navy", ha="center", va="center",
             transform=ax2.transAxes,fontsize=30, fontweight='bold',
             rotation = 45, alpha = 0.3, bbox = dict(facecolor='none', edgecolor='none'))

ax2.set_xlim(bin_min, bin_max)
ax2.set_xlabel("E/p", fontsize=12)
ax2.set_ylabel("Counts", fontsize=12)
ax2.set_title(f"Fitted E/p", fontsize=14)
ax2.grid(alpha=0.5)
ax2.legend()

# --------------------------------------------------------------------------
# Save the combined figure
# --------------------------------------------------------------------------
plt.savefig(f"{file_type}/{file_type}_run_{runnum}_pcal.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print(f"{runnum}\t{mean_fit}\t{mean_err}\t{sigma_fit}\t{sigma_err}\t"
      f"{fit_bin_min}\t{fit_bin_max}\t{fit_bin_avg}\t{fit_bin_sum}")

