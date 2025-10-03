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

runnum = input(f"Input the run number you wish to analyze:")

# Directory information, modify as needed
root_directory = f"/volatile/hallc/c-rsidis/cmorean/replay_pass0a/ROOTfiles"
coin_pattern = f"coin_replay_production_{runnum}_-1.root"
shms_pattern = f"hms_coin_replay_production_{runnum}_-1.root"

d_calo_fp = 292.64 # distance from focal plane to calorimeter face

branches = [
    "P.dc.x_fp",
    "P.dc.y_fp",
    "P.dc.xp_fp",
    "P.dc.yp_fp",
    "P.gtr.dp",
    "P.ngcer.npeSum",
    "P.cal.etottracknorm"
]


try:
    input_root_filepath = f"{root_directory}/{coin_pattern}"
    df = pd.DataFrame(uproot.open(input_root_filepath)["T"].arrays(branches, library = "np"))
    file_type = "COIN"
except OSError:
    try:
        input_root_filepath = f"{root_directory}/{shms_pattern}"
        df = pd.DataFrame(uproot.open(input_root_filepath)["T"].arrays(branches, library = "np"))
        file_type = "SHMS"
    except OSError:
        print(f"Run {runnum} not valid for pcal calibration; skipping...")
        sys.exit(1)
        
print(f"Found {file_type} run {runnum}.  Dataframe {len(df)} rows.  Processing...")

cuts = ((df["P.gtr.dp"] > -15) & (df["P.gtr.dp"] < 27) & (df["P.ngcer.npeSum"] > 2))

df_cut = df[cuts].copy()

xcalo = ((df_cut["P.dc.x_fp"]) + (df_cut["P.dc.xp_fp"])*d_calo_fp)

ycalo = ((df_cut["P.dc.y_fp"]) + (df_cut["P.dc.yp_fp"])*d_calo_fp)

weight = df_cut["P.cal.etottracknorm"]

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
# Plotting figures, disabling most of them but will leave in code for debugging
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
# Plotting
# -----------------------------------------------------------------------------

# Setting up custom color mapping
data = hCaloPosNormU.view()
vmin, vmax = np.nanmin(data), np.nanmax(data)
norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)

mid = norm(1.0)

cmap = mplcolors.LinearSegmentedColormap.from_list(
    "custom_cmap",
    [(0.0, "white"), (mid, "navy"), (0.5, "darkorange"), (1.0, "red")]
)

# Plotting normalized e/p per track
plt.figure(figsize=(8, 6))
plt.imshow(
    hCaloPosNormU.view().T,
    origin="lower",
    extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
    aspect="auto",
    cmap = cmap
)
fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    hCaloPosNormU.view().T,
    origin="lower",
    extent=[xrange[0], xrange[1], yrange[0], yrange[1]],
    aspect="auto",
    cmap= cmap
)

ax.set_title(f"Run {runnum} Normalized E/p per Track at SHMS Calorimeter")

imvmin, imvmax = im.get_clim()
cbar = fig.colorbar(im, ax=ax, pad = 0.01, format='%.1f')
cbar.ax.tick_params(labelsize=8)  # Consistent font size and spacing


# Plotting the grid
numrows, numcols = 16, 14
blockspacing = 9

startrow = -blockspacing * numrows / 2.0
startcol = -blockspacing * numcols / 2.0

for i in range(numrows + 1):
    y = startrow + i * blockspacing
    ax.hlines(y, startcol, startcol + blockspacing * numcols, colors='silver', linewidth=1.0, alpha = 0.6)

for i in range(numcols + 1):
    x = startcol + i * blockspacing
    ax.vlines(x, startrow, startrow + blockspacing * numrows, colors='silver', linewidth=1.0, alpha = 0.6)

plt.savefig(f"{file_type}/{file_type}_run_{runnum}_etottracknorm_at_cal.png", bbox_inches="tight", dpi=300)
plt.close()

# Plotting e/p

bin_min, bin_max, bin_num = 0, 2, 200
data_bins = np.linspace(bin_min, bin_max, bin_num + 1)

counts, bin_edges = np.histogram(df_cut["P.cal.etottracknorm"], bins=data_bins)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Fitting range
fit_min, fit_max = 0.9, 1.1
fit_mask = (bin_centers >= fit_min) & (bin_centers <= fit_max)

x_fit = bin_centers[fit_mask]
y_fit = counts[fit_mask]

# Mask out zero bins because log-likelihood or fit weighting can fail
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
popt, _ = curve_fit(gaussian, x_fit, y_fit, p0=[amp_guess, mean_guess, sigma_guess])

plt.figure()
plt.hist(df_cut["P.cal.etottracknorm"], bins=data_bins, histtype='step', color='red', label='e/p')
plt.xlim(bin_min, bin_max)
plt.xlabel("E/p", fontsize=12)
plt.ylabel("Counts", fontsize=12)
plt.title(f"Run {runnum} Fitted E/p", fontsize=16)
plt.grid(alpha=0.5)

# Plot the gaussian fit
x_plot = np.linspace(fit_min, fit_max, 500)
plt.plot(x_plot, gaussian(x_plot, *popt), 'b-', label='Gaussian fit')

plt.legend()
plt.tight_layout()
plt.savefig("Testing.png")
