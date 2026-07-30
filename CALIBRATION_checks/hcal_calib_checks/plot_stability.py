#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter  #Used for tick mark additions
import mplhep

fit_results_filepath = "CSVs/FIT_hcal_results.csv"
outdir = "STABILITY_plots"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(fit_results_filepath, comment="#")
df.columns = [c.strip() for c in df.columns]

for col in ["runnum", "hms_p", "hms_th", "fit_mean", "mean_err", "fit_sigma", "sigma_err", "bin_min", "bin_max", "bin_avg", "bin_total"]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

df["run_type"] = df["run_type"].astype(str).str.strip()
df = df.dropna(subset=["runnum", "run_type", "hms_p", "fit_mean", "mean_err"]).copy()

df["runnum_int"] = df["runnum"].astype(int)
df["momentum"] = df["hms_p"].round(3)
df["residual"] = df["fit_mean"] - 1.0
df = df.sort_values(["momentum", "runnum_int"])

available_momenta = sorted(df["momentum"].unique())
available_types = sorted(df["run_type"].unique())

# mom_colors = {mom: plt.cm.Set1.colors[i % 10] for i, mom in enumerate(available_momenta)}
# cmap = plt.cm.gnuplot
# colors = cmap(np.linspace(0.1, 0.9, len(available_momenta)))

# mom_colors = {
#     mom: colors[i]
#     for i, mom in enumerate(available_momenta)
# }

angle_values = sorted(df["hms_th"].dropna().unique()) if "hms_th" in df.columns else []
# angle_colors = {ang: plt.cm.Set1.colors[i % 20] for i, ang in enumerate(angle_values)}

def plot_all_momenta(plot_df, outfile):
    plt.style.use(mplhep.style.ROOT)
    plt.rcParams.update({"figure.titlesize": 14,
                         "axes.titlesize": 14,
                         "axes.labelsize": 14,
                         "legend.fontsize": 12,
                         "xtick.labelsize": 12,
                         "ytick.labelsize": 12})
    
    fig, (ax_top, ax_bot) = plt.subplots(nrows=2, ncols=1, figsize=(8, 6), gridspec_kw={"height_ratios": [3, 1]}, constrained_layout = True)

    for mom in sorted(plot_df["momentum"].unique()):
        subdf = plot_df[plot_df["momentum"] == mom].sort_values("runnum_int")
        ax_top.errorbar(subdf["runnum_int"].to_numpy(), subdf["fit_mean"].to_numpy(), yerr=subdf["mean_err"].to_numpy(), fmt="o", capsize=0, elinewidth=0.5, markersize=3, label=f"{float(mom):.3f} GeV")

    ax_top.axhline(1.0, color="black", linestyle="-", linewidth=1.5)

    subdf = plot_df.sort_values("runnum_int")

    x = subdf["runnum_int"].to_numpy()
    y = subdf["fit_mean"].to_numpy()
    s = subdf["fit_sigma"].to_numpy()

    ax_top.fill_between(x, y - 3.0 * s, y + 3.0 * s, alpha=0.36, zorder=0, interpolate=False, label = "$\pm 3\sigma$ band")
    ax_top.set_ylabel(r"H.cal.ettotracknorm Fit")
    ax_top.set_xlabel("Run Number")
    ax_top.grid(True)
    ax_top.legend(loc="upper right", frameon = True, fancybox = True, framealpha = 0.6, edgecolor = "gray")

    res = plot_df["residual"].dropna()
    ax_bot.hist(res, bins=100, histtype="stepfilled", alpha=0.36)
    ax_bot.axvline(0, color="black", linestyle="-", linewidth=1.5)

    sigma = res.std()
    if np.isfinite(sigma) and sigma > 0:
        for n in [1, 2, 3]:
            ax_bot.axvline(n * sigma, color="red", linestyle="--", alpha=0.6)
            ax_bot.axvline(-n * sigma, color="red", linestyle="--", alpha=0.6)

    ax_bot.set_xlabel("Residual (Fit Mean − 1.0)")
    ax_bot.set_ylabel("Count")
    ax_bot.grid(True)

    ax_top.set_ylim(0.7,)

    fig.suptitle(f"HMS Calorimeter Stability (All Settings)\nRun Types: {', '.join(available_types)}", fontsize=15)
    # fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)

def plot_single_momentum(mom, plot_df, outfile):
    plt.style.use(mplhep.style.ROOT)
    plt.rcParams.update({"figure.titlesize": 14,
                         "axes.titlesize": 14,
                         "axes.labelsize": 14,
                         "legend.fontsize": 12,
                         "xtick.labelsize": 12,
                         "ytick.labelsize": 12})
    
    fig, (ax_top, ax_bot) = plt.subplots(nrows=2, ncols=1, figsize=(8, 6), gridspec_kw={"height_ratios": [3, 1]},constrained_layout=True)

    for ang in sorted(plot_df["hms_th"].dropna().unique()):
        subdf = plot_df[plot_df["hms_th"] == ang].sort_values("runnum_int")
        ax_top.errorbar(subdf["runnum_int"].to_numpy(), subdf["fit_mean"].to_numpy(), yerr=subdf["mean_err"].to_numpy(), fmt="o", capsize=0, elinewidth=0.5, markersize=3, label=f"{float(ang):.3f} deg")

    ax_top.axhline(1.0, color="black", linestyle="-", linewidth=1.5)
    ax_top.set_ylabel(r"H.cal.etottracknorm Fit")
    ax_top.set_xlabel("Run Number")
    ax_top.grid(True)
    ax_top.legend(loc="upper right", frameon = True, fancybox = True, framealpha = 0.6, edgecolor = "gray")

    for ang in sorted(plot_df["hms_th"].dropna().unique()):
        subres = plot_df.loc[plot_df["hms_th"] == ang, "residual"].dropna()
        ax_bot.hist(subres, bins=100, histtype="stepfilled", alpha = 0.6, linewidth=1.2, label=f"{float(ang):.3f} deg")

    ax_bot.axvline(0, color="black", linestyle="-", linewidth=1.5)
    ax_bot.set_xlabel("Residual (Fit Mean − 1.0)")
    ax_bot.set_ylabel("Count")
    ax_bot.grid(True)
    ax_bot.legend(loc="upper right", frameon = True, fancybox = True, framealpha = 0.6, edgecolor = "gray")

    fig.suptitle(f"HMS Calorimeter Stability ({float(mom):.3f} GeV)\nRun Types: {', '.join(sorted(plot_df['run_type'].unique()))}")
    fig.savefig(outfile, dpi=300)
    plt.close(fig)

plot_all_momenta(df, os.path.join(outdir, "all_momenta_stability_hcal.png"))

for mom in available_momenta:
    momdf = df[df["momentum"] == mom].copy()
    plot_single_momentum(mom, momdf, os.path.join(outdir, f"stability_hcal_{float(mom):.3f}GeV.png"))

print(f"Wrote {1 + len(available_momenta)} plots to {outdir}")
