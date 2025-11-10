#!/usr/bin/env python3

import sys, os
import matplotlib
matplotlib.use('Agg')
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import boost_histogram as bh
import matplotlib.colors as mplcolors
from scipy.optimize import curve_fit

# --------------------------------------------------------------------------
# User inputs, directories, branches for plotting scripts; modify as needed
# --------------------------------------------------------------------------
if len(sys.argv) > 1:
    runnum = sys.argv[1]
else:
    runnum = input("Input the run number you wish to analyze: ")

root_directory = f"/volatile/hallc/c-rsidis/relder/ROOTfiles"
coin_pattern = f"coin_replay_production_{runnum}_-1.root"
shms_pattern = f"shms_coin_replay_production_{runnum}_-1.root"
hms_pattern = f"hms_coin_replay_production_{runnum}_-1.root"
auxfiles_runlist_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat"

d_calo_fp = 292.64 # distance from focal plane to calorimeter face

branches = ["P.dc.x_fp", "P.dc.y_fp", "P.dc.xp_fp", "P.dc.yp_fp",
            "P.gtr.dp", "P.ngcer.npeSum", "P.cal.etottracknorm"]

# --------------------------------------------------------------------------
# Trying to open root files, basically filters out HMS replays
# --------------------------------------------------------------------------
try:
    input_root_filepath = f"{root_directory}/{coin_pattern}"
    if not os.path.exists(input_root_filepath):
        raise FileNotFoundError
    df = uproot.open(input_root_filepath)["T"].arrays(branches, library="pd")
    file_type = "COIN"

except FileNotFoundError:
    try:
        input_root_filepath = f"{root_directory}/{shms_pattern}"
        if not os.path.exists(input_root_filepath):
            raise FileNotFoundError
        df = uproot.open(input_root_filepath)["T"].arrays(branches, library="pd")
        file_type = "SHMS"

    except FileNotFoundError:
        input_root_filepath = f"{root_directory}/{hms_pattern}"
        if os.path.exists(input_root_filepath):
            print(f"WARNING:\tRun {runnum} is an HMS replay; skipping...")
        else:
            print(f"ERROR:\tRun {runnum} is missing a ROOTfile; skipping...")
        sys.exit(1)

    except uproot.exceptions.KeyInFileError as e:
        print(f"ERROR:\tRun {runnum} is an SHMS file missing branch '{e.key}'; skipping...")
        sys.exit(1)

except uproot.exceptions.KeyInFileError as e:
    print(f"ERROR:\tRun {runnum} is a COIN file missing branch '{e.key}'; skipping...")
    sys.exit(1)

# --------------------------------------------------------------------------
# Setting up variables for plotting
# --------------------------------------------------------------------------
cuts = ((df["P.gtr.dp"] > -15) & (df["P.gtr.dp"] < 27) & (df["P.ngcer.npeSum"] > 2)) & (df["P.cal.etottracknorm"] > 0)

df_cut = df[cuts]

xcalo = ((df_cut["P.dc.x_fp"]) + (df_cut["P.dc.xp_fp"])*d_calo_fp)

ycalo = ((df_cut["P.dc.y_fp"]) + (df_cut["P.dc.yp_fp"])*d_calo_fp)

weight = df_cut["P.cal.etottracknorm"]

finite_mask = np.isfinite(xcalo) & np.isfinite(ycalo) & np.isfinite(weight)

xcalo = xcalo[finite_mask]

ycalo = ycalo[finite_mask]

weight = weight[finite_mask]

xbins, ybins = 200, 200

xmin, xmax = -60, 60

ymin, ymax = -58, 58

xrange = (xmin - 10, xmax + 10)

yrange = (ymin - 10, ymax + 10)

hCaloPos = bh.Histogram(bh.axis.Regular(xbins, xrange[0], xrange[1]), bh.axis.Regular(ybins, yrange[0], yrange[1]))
hCaloPos.fill(ycalo, xcalo)

hCaloPosWtU = bh.Histogram(bh.axis.Regular(xbins, xrange[0], xrange[1]), bh.axis.Regular(ybins, yrange[0], yrange[1]))
hCaloPosWtU.fill(ycalo,xcalo,weight=weight)

counts_weighted = hCaloPosWtU.view()
counts_unweighted = hCaloPos.view()

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
data = hCaloPosNormU.view()
vmin = 0
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

ax1.set_title(f"Normalized E/p per Track at Calorimeter", fontsize=14)

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
bin_min, bin_max, bin_num = 0, 2, 200
data_bins = np.linspace(bin_min, bin_max, bin_num + 1)
counts, bin_edges = np.histogram(df_cut["P.cal.etottracknorm"], bins=data_bins)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

fit_min, fit_max = 0.90, 1.10
fit_mask = (bin_centers >= fit_min) & (bin_centers <= fit_max)

x_fit = bin_centers[fit_mask]
y_fit = counts[fit_mask]

nonzero = y_fit > 0
x_fit = x_fit[nonzero]
y_fit = y_fit[nonzero]

# Guess for initial fit ONLY within fit range
if y_fit.size == 0:
    print(f"WARNING:\tRun {runnum} contains no events in range; skipping...")
    sys.exit(2)
# if np.sum(y_fit) < 50:
#     print(f"WARNING:\tRun {runnum} does not contain enough events in range; skipping...")
#     sys.exit(2)
amp_guess = np.max(y_fit)
mean_guess = np.sum(x_fit * y_fit) / np.sum(y_fit)
sigma_guess = np.sqrt(np.sum(y_fit * (x_fit - mean_guess) ** 2) / np.sum(y_fit))

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

# Setting a flag for failures in the fit.  If unable to find a fit, we'll still make the plot, just with printed text on top.
fit_failed = False
try:
    popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=[amp_guess, mean_guess, sigma_guess])
    amp_fit, mean_fit, sigma_fit = popt
    amp_err, mean_err, sigma_err = np.sqrt(np.diag(pcov))
except:
    fit_failed = True
    print(f"WARNING:\tRun {runnum} was unable to find a Gaussian fit; skipping...")
    amp_fit = mean_fit = sigma_fit = amp_err = mean_err = sigma_err = np.nan

# Now plotting to the second pad
ax2.hist(df_cut["P.cal.etottracknorm"], bins=data_bins, histtype='step', color='red', label='E/p')
x_plot = np.linspace(fit_min, fit_max, 500)

if not fit_failed:
    x_plot = np.linspace(fit_min, fit_max, 500)
    ax2.plot(x_plot, gaussian(x_plot, *popt), color='navy', label='Gaussian fit')
    ax2.text(0.05, 0.95,
             f"Mean = {mean_fit:.4f}\nSigma = {sigma_fit:.4f}",
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
plt.close('all')

print(f"{runnum}\t{mean_fit}\t{mean_err}\t{sigma_fit}\t{sigma_err}") # Don't uncomment this line!  This helps for PLOT_... script.

