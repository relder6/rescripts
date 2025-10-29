#!/usr/bin/env python3

import sys
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

root_directory = f"/work/hallc/c-rsidis/skimfiles/pass0"
coin_pattern = f"skimmed_coin_replay_production_{runnum}_-1.root"
hms_pattern = f"skimmed_hms_coin_replay_production_{runnum}_-1.root"
auxfiles_runlist_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat"

d_calo_fp = 338.69 # distance from focal plane to calorimeter face

branches = ["H_dc_x_fp", "H_dc_y_fp", "H_dc_xp_fp", "H_dc_yp_fp",
            "H_gtr_dp", "H_cer_npeSum", "H_cal_etottracknorm"]

# --------------------------------------------------------------------------
# Trying to open root files, basically filters out SHMS replays
# --------------------------------------------------------------------------
try:
    input_root_filepath = f"{root_directory}/{coin_pattern}"
    df = pd.DataFrame(uproot.open(input_root_filepath)["T"].arrays(branches, library = "np"))
    file_type = "COIN"
except OSError:
    try:
        input_root_filepath = f"{root_directory}/{hms_pattern}"
        df = pd.DataFrame(uproot.open(input_root_filepath)["T"].arrays(branches, library = "np"))
        file_type = "HMS"
    except OSError:
        print(f"Run {runnum} not valid for hcal calibration; skipping...")
        sys.exit(1)
        
print(f"Found {file_type} run {runnum}.  Dataframe {len(df)} rows.  Processing...")

# --------------------------------------------------------------------------
# Setting up variables for plotting
# --------------------------------------------------------------------------
cuts = ((df["H_gtr_dp"].between(-8,8)) & (df["H_cer_npeSum"] > 1.5))

df_cut = df[cuts].copy()

xcalo = ((df_cut["H_dc_x_fp"]) + (df_cut["H_dc_xp_fp"])*d_calo_fp)

ycalo = ((df_cut["H_dc_y_fp"]) + (df_cut["H_dc_yp_fp"])*d_calo_fp)

weight = df_cut["H_cal_etottracknorm"]

finite_mask = np.isfinite(xcalo) & np.isfinite(ycalo) & np.isfinite(weight)

xcalo = xcalo[finite_mask]

ycalo = ycalo[finite_mask]

weight = weight[finite_mask]

xbins, ybins = 100, 100

xmin, xmax = -65.4, 54.6

ymin, ymax = -30, 30

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
# Plotting figures, disabling most of them but will leave in code for debugging
# -----------------------------------------------------------------------------
# Plotting figures, disabling most of them but will leave in code for debugging
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

# # Plot and save hExitPos
# plt.figure(figsize=(8, 6))
# plt.imshow(
#     hExitPos.view().T,
#     origin="lower",
#     extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
#     aspect="auto"
# )
# plt.xlabel("xexit")
# plt.ylabel("yexit")
# plt.title("2D Histogram: yexit vs xexit")
# plt.colorbar(label="Counts")
# plt.savefig("hExitPos.png", bbox_inches="tight", dpi=150)
# plt.close()

# # Plotting and saving hCaloPosWtU
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

# # Plotting and saving hCaloPosNormU
# plt.figure(figsize=(8, 6))
# plt.imshow(
#     hCaloPosNormU.view().T,
#     origin="lower",
#     extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
#     aspect="auto",
#     cmap = "afmhot_r"
# )
# plt.xlabel("xnorm")
# plt.ylabel("ynorm")
# plt.title("ynorm v xnorm")
# plt.colorbar(label="Counts")
# plt.savefig("hCaloPosNormU.png", bbox_inches="tight", dpi=150)
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
vmax = vmax_candidate if vmax_candidate > vmin else vmin + 1  # avoid zero range here

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

# Now adding the block grid overlay,
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
counts, bin_edges = np.histogram(df_cut["H_cal_etottracknorm"], bins=data_bins)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

fit_min, fit_max = 0.9, 1.1
fit_mask = (bin_centers >= fit_min) & (bin_centers <= fit_max)

x_fit = bin_centers[fit_mask]
y_fit = counts[fit_mask]

nonzero = y_fit > 0
x_fit = x_fit[nonzero]
y_fit = y_fit[nonzero]

# Guess for initial fit ONLY within fit range
amp_guess = np.max(y_fit)
mean_guess = np.sum(x_fit * y_fit) / np.sum(y_fit)
sigma_guess = np.sqrt(np.sum(y_fit * (x_fit - mean_guess) ** 2) / np.sum(y_fit))

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

from scipy.optimize import curve_fit
popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=[amp_guess, mean_guess, sigma_guess])
amp_fit, mean_fit, sigma_fit = popt
amp_err, mean_err, sigma_err = np.sqrt(np.diag(pcov))

# Extracting here the fit results
amp_fit, mean_fit, sigma_fit = popt

# Now plotting to the second pad
ax2.hist(df_cut["H_cal_etottracknorm"], bins=data_bins, histtype='step', color='red', label='E/p')
x_plot = np.linspace(fit_min, fit_max, 500)
ax2.plot(x_plot, gaussian(x_plot, *popt), color = 'navy', label='Gaussian fit')

ax2.set_xlim(bin_min, bin_max)
ax2.set_xlabel("E/p", fontsize=12)
ax2.set_ylabel("Counts", fontsize=12)
ax2.set_title(f"Fitted E/p", fontsize=14)
ax2.grid(alpha=0.5)
ax2.legend()

ax2.text(
    0.05, 0.95,
    f"Mean = {mean_fit:.4f}\nSigma = {sigma_fit:.4f}",
    transform=ax2.transAxes,
    verticalalignment='top',
    fontsize=10,
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='black')
)

# --------------------------------------------------------------------------
# Save the combined figure
# --------------------------------------------------------------------------

plt.savefig(f"{file_type}/{file_type}_run_{runnum}_hcal.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"{runnum}\t{mean_fit}\t{mean_err}\t{sigma_fit}\t{sigma_err}") # Don't uncomment this line!  This helps for PLOT_... script.

