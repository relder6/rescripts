#!/usr/bin/env python3

import sys, os, argparse, uproot, csv
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import boost_histogram as bh
import matplotlib.colors as mplcolors
from scipy.optimize import curve_fit

# --------------------------------------------------------------------------
# Check file locations before running,
# --------------------------------------------------------------------------
root_directory = f"/w/hallc-scshelf2102/c-rsidis/skimfiles/pass0p1"
bigtable_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_bigtable_pass0p1.csv"
outfile = "CSVs/FIT_pcal_ep_results.csv"

# Parsing optional command line input for selected run numbers
parser = argparse.ArgumentParser()
parser.add_argument("runs", nargs="*", type=int, help="Optional run numbers")
args = parser.parse_args()
selected_runs = set(args.runs) if args.runs else None

data = np.genfromtxt(bigtable_filepath,delimiter=",",names=True,dtype=None,encoding=None)

type_mask = np.isin(data["run_type"], ["SHMSDIS", "PI-SIDIS"])
polarity_mask = data["shms_p"] < 0
mask = type_mask & polarity_mask

runnums = data["run"][mask]
run_types = data["run_type"][mask]
shms_ps = data["shms_p"][mask]

if selected_runs is not None:
    run_mask = np.isin(runnums, list(selected_runs))
    runnums = runnums[run_mask]
    run_types = run_types[run_mask]
    shms_ps = shms_ps[run_mask]

print(f"Found {len(runnums)} runs")

# Defining some common things that will be used in every run analysis,
d_calo_fp = 292.64 # distance from focal plane to calorimeter face

branches = ["P_dc_x_fp", "P_dc_y_fp", "P_dc_xp_fp", "P_dc_yp_fp", "P_gtr_dp", "P_hgcer_npeSum", "P_ngcer_npeSum", "P_cal_etottracknorm", "P_gtr_beta"]
data_cut = ("(P_gtr_dp > -15) & (P_gtr_dp < 27) & (P_hgcer_npeSum > 1.5) & (P_ngcer_npeSum > 1.5) & (P_gtr_beta > 0.8) & (P_gtr_beta < 1.2) & (P_cal_etottracknorm > 0)")

xbins, ybins = 100, 100
xmin, xmax = -60, 60
ymin, ymax = -58, 58
xrange = (xmin - 10, xmax + 10)
yrange = (ymin - 10, ymax + 10)

bin_min, bin_max, bin_num = 0, 2, 100
data_bins = np.linspace(bin_min, bin_max, bin_num + 1)

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

def safe_norm_pos(val, vmin, vmax, prev_pos):
    if vmax == vmin:
        return 1.0
    pos = (val - vmin) / (vmax - vmin)
    pos = np.clip(pos, 0.0, 1.0)
    return max(pos, prev_pos)

grid_numrows, grid_numcols, blockspacing = 16, 14, 9
startrow, startcol = -blockspacing * grid_numrows / 2.0, -blockspacing * grid_numcols / 2.0
grid_ys = [startrow + i * blockspacing for i in range(grid_numrows + 1)]
grid_xs = [startcol + i * blockspacing for i in range(grid_numcols + 1)]

with open(outfile, "w", newline="") as csvfile:
    flush_every = 10
    written = 0
    writer = csv.writer(csvfile)

    writer.writerow(["runnum","run_type","shms_p","fit_mean","mean_err",
                     "fit_sigma","sigma_err","bin_min","bin_max",
                     "bin_avg","bin_total"])

    for runnum, run_type, shms_p in zip(runnums, run_types, shms_ps):
        print(f"Processing run {runnum} ({run_type})")
        if run_type == "SHMSDIS":
            skimfile = f"{root_directory}/skimmed_shms_coin_replay_production_{runnum}_-1.root"
        elif run_type == "PI-SIDIS":
            skimfile = f"{root_directory}/skimmed_coin_replay_production_{runnum}_-1.root"

        if not os.path.exists(skimfile):
            print(f"ERROR:\tRun {runnum} is missing a ROOTfile; skipping...")
            continue

        try:
            with uproot.open(skimfile) as rootfile:
                tree = rootfile["T"]
                arrays = tree.arrays(branches, cut=data_cut, library="np")
        except uproot.exceptions.KeyInFileError as e:
            print(f"ERROR:\tRun {runnum} is a {run_type} file missing branch '{e.key}'; skipping...")
            sys.exit(1) 
    
        # --------------------------------------------------------------------------
        # Setting up variables for plotting
        # --------------------------------------------------------------------------
        ep = arrays["P_cal_etottracknorm"]
        
        xcalo = (arrays["P_dc_x_fp"]) + (arrays["P_dc_xp_fp"])*d_calo_fp

        ycalo = (arrays["P_dc_y_fp"]) + (arrays["P_dc_yp_fp"])*d_calo_fp
        
        weight = ep
        
        finite_mask = np.isfinite(xcalo) & np.isfinite(ycalo) & np.isfinite(weight)
        
        xcalo = xcalo[finite_mask]
        
        ycalo = ycalo[finite_mask]
        
        weight = weight[finite_mask]

        hCaloPos = bh.Histogram(bh.axis.Regular(xbins, xrange[0], xrange[1]),
                                bh.axis.Regular(ybins, yrange[0], yrange[1]))

        hCaloPosWt = bh.Histogram(bh.axis.Regular(xbins, xrange[0], xrange[1]),
                                  bh.axis.Regular(ybins, yrange[0], yrange[1]),
                                  storage = bh.storage.Weight())

        hCaloPos.fill(ycalo, xcalo)
        hCaloPosWt.fill(ycalo, xcalo, weight=weight)

        counts_unweighted = hCaloPos.values()
        counts_weighted = hCaloPosWt.values()

        total_counts = counts_unweighted.sum()
        if total_counts == 0:
            print(f"WARNING:\tRun {runnum} has zero tracks; skipping...")
            continue

        with np.errstate(divide='ignore', invalid='ignore'):
            hCaloPosNormU = np.divide(counts_weighted, counts_unweighted)
            hCaloPosNormU[~np.isfinite(hCaloPosNormU)] = 0  # set inf/NaN to 0

        # -----------------------------------------------------------------------------
        # Getting the figure ready to place both plots
        # -----------------------------------------------------------------------------
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (14, 6), constrained_layout = True)
        fig.suptitle(f"{run_type} Run {runnum}", fontsize=16, fontweight="bold")

        # -----------------------------------------------------------------------------
        # Plotting normalized e/p per track at calorimeter
        # -----------------------------------------------------------------------------
        data = hCaloPosNormU
        vmin = 0.0
        vmax_candidate = np.nanmax(data)
        vmax = vmax_candidate if vmax_candidate > vmin else vmin + 1  # avoiding zero range here

        pos_0 = 0.0
        pos_1 = safe_norm_pos(1, vmin, vmax, pos_0)
        pos_2 = safe_norm_pos(2, vmin, vmax, pos_1)
        pos_max = 1.0

        cmap = mplcolors.LinearSegmentedColormap.from_list(
            "custom_cmap",[(pos_0, "white"),
                           (pos_1, "navy"),
                           (pos_2, "orange"),
                           (pos_max, "red")])

        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)

        im = ax1.imshow(data.T,origin="lower",extent=[xrange[0], xrange[1], yrange[0], yrange[1]],aspect="auto",cmap=cmap,norm=norm)
        
        cbar = fig.colorbar(im, ax=ax1, format='%.1f', pad=0.01)

        ticks = np.arange(0, int(np.floor(vmax)) + 0.5, 1)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f"{tick:.1f}" for tick in ticks])
        cbar.ax.tick_params(labelsize=8)

        ax1.set_title(f"Normalized E/p per Track at Calorimeter", fontsize = 14)

        for y in grid_ys:
            ax1.hlines(y, startcol, startcol + blockspacing * grid_numcols, colors="silver", linewidth=1.0, alpha=0.6)
        for x in grid_xs:
            ax1.vlines(x, startrow, startrow + blockspacing * grid_numrows, colors="silver", linewidth=1.0, alpha=0.6)
        # -----------------------------------------------------------------------------
        # Plotting fitted distribution of e/p
        # -----------------------------------------------------------------------------
        counts, bin_edges = np.histogram(arrays["P_cal_etottracknorm"], bins=data_bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Initializing fit_failed flag, several checks will be made against this to avoid wasting time fitting nothing
        fit_failed = False

        # Need a broad search window here; want to capture peaks that are poorly fit (--> bad calibration!)
        search_min, search_max = 0.80, 1.20
        min_peak_height = 15
        search_mask = (bin_centers >= search_min) & (bin_centers <= search_max)

        x_search = bin_centers[search_mask]
        y_search = counts[search_mask]

        nonzero = y_search > 0

        x_search = x_search[nonzero]
        y_search = y_search[nonzero]

        fit_bin_min = fit_bin_max = fit_bin_avg = fit_bin_sum = np.nan
        mean_fit = mean_err = sigma_fit = sigma_err = np.nan

        if x_search.size == 0:
            print(f"WARNING:\tRun {runnum} no entries in search window; skipping...")
            fit_failed = True

        # Here I'm finding the peak in the broad search window, and only keeping 0.15 around that,

        else:   
            peak_idx = np.argmax(y_search)
            peak_x = x_search[peak_idx]
            peak_y = y_search[peak_idx]

            if peak_y < min_peak_height:
                fit_failed = True
                print(f"WARNING:\tRun {runnum} peak too small ({peak_y}); skipping...")
            else:
                fit_width = 0.15
                local_mask = np.abs(x_search - peak_x) <= fit_width
                x_fit = x_search[local_mask]
                y_fit = y_search[local_mask]

                if x_fit.size == 0:
                    fit_failed = True
                    print(f"WARNING:\tRun {runnum} contains no events in range; skipping...")
                else:
                    keep_above = 0.30
                    peak_height = np.max(y_fit)
                    top_mask = y_fit >= keep_above * peak_height
                    x_fit = x_fit[top_mask]
                    y_fit = y_fit[top_mask]

                    if x_fit.size < 3:
                        fit_failed = True
                        print(f"WARNING:\tRun {runnum} not enough bins around peak; skipping...")

                    else:
                        fit_bin_min = np.min(x_fit)
                        fit_bin_max = np.max(x_fit)
                        fit_bin_avg = np.mean(x_fit)
                        fit_bin_sum = np.sum(y_fit)
    
                        amp_guess = np.max(y_fit)
                        mean_guess = np.sum(x_fit * y_fit) / np.sum(y_fit)
                        sigma_guess = np.sqrt(np.sum(y_fit * (x_fit - mean_guess) ** 2) / np.sum(y_fit))

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
        ax2.hist(arrays["P_cal_etottracknorm"], bins=data_bins, histtype='step', color='red', label='E/p')

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
        fig.savefig(f"{run_type}/{run_type}_run_{runnum}_pcal.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        writer.writerow([runnum,run_type,shms_p,mean_fit,mean_err,sigma_fit,sigma_err,fit_bin_min,fit_bin_max,fit_bin_avg,fit_bin_sum])
        written += 1
        if written % flush_every == 0:
            csvfile.flush()
    csvfile.flush()
            
