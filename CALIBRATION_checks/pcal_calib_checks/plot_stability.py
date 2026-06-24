#!/usr/bin/env python3

import os, re, sys, csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fit_results_filepath = "CSVs/FIT_pcal_ep_results.csv"
outdir = "STABILITY_plots"
os.makedirs(outdir, exist_ok=True)

rows = []

df = pd.read_csv(fit_results_filepath, comment="#")
df.columns = [c.strip() for c in df.columns]

for col in ["runnum", "shms_p", "fit_mean", "mean_err", "fit_sigma", "sigma_err", "bin_min", "bin_max", "bin_avg", "bin_total"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df["run_type"] = df["run_type"].astype(str).str.strip()
df = df.dropna(subset=["runnum", "run_type", "shms_p", "fit_mean", "mean_err"]).copy()

df["runnum_int"] = df["runnum"].astype(int)
df["momentum"] = df["shms_p"].round(3)
df = df.sort_values(["momentum", "runnum_int"])
df["residual"] = df["fit_mean"] - 1.0

available_momenta = sorted(df["momentum"].unique())
available_types = sorted(df["run_type"].unique())

type_markers = {"SHMSDIS": "o", "PI-SIDIS": "s"}
default_markers = ["o", "s", "D", "^", "v", "x", "P", "*"]
for i, rt in enumerate(available_types):
    type_markers.setdefault(rt, default_markers[i % len(default_markers)])

colors = plt.cm.tab10.colors
momentum_colors = {mom: colors[i % len(colors)] for i, mom in enumerate(available_momenta)}

def make_plot(plot_df, title_text, outfile):
    fig, (ax_top, ax_bot) = plt.subplots(
        nrows=2, ncols=1, figsize=(8, 8),
        gridspec_kw={"height_ratios": [3, 1]}
    )

    for mom in sorted(plot_df["momentum"].unique()):
        momdf = plot_df[plot_df["momentum"] == mom]
        color = momentum_colors[mom]
        for rt in sorted(momdf["run_type"].unique()):
            subdf = momdf[momdf["run_type"] == rt].sort_values("runnum_int")
            ax_top.errorbar(
                subdf["runnum_int"].to_numpy(),
                subdf["fit_mean"].to_numpy(),
                yerr=subdf["mean_err"].to_numpy(),
                fmt=type_markers.get(rt, "o"),
                color=color,
                capsize=0,
                elinewidth=0.5,
                markersize=3,
                label=f"{mom:g} GeV, {rt}"
            )

    ax_top.axhline(1.0, color="navy", linestyle="--", linewidth=1.2, label="y = 1")
    ax_top.set_ylabel(r"$\mu_{\rm fit}$, P.cal.etottracknorm", fontsize=13)
    ax_top.set_xlabel("Run Number", fontsize=13)
    ax_top.legend(title="Momentum, Run Type", fontsize=9, ncol=2)
    ax_top.grid(True, linestyle="--", alpha=0.6)

    res = plot_df["residual"].dropna()
    ax_bot.hist(res, bins=100, histtype="stepfilled", alpha=0.7, edgecolor="black")
    ax_bot.axvline(0, color="black", linestyle="--", linewidth=1)

    sigma = res.std()
    if np.isfinite(sigma) and sigma > 0:
        for n in [1, 2, 3]:
            ax_bot.axvline(n * sigma, color="red", linestyle="--", alpha=0.6)
            ax_bot.axvline(-n * sigma, color="red", linestyle="--", alpha=0.6)

    ax_bot.set_xlabel("Residual (Fit Mean − 1.0)", fontsize=13)
    ax_bot.set_ylabel("Count", fontsize=13)
    ax_bot.grid(True, linestyle="--", alpha=0.6)

    fig.suptitle(title_text, fontsize=15)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)

make_plot(
    df,
    "SHMS Calorimeter Stability (All Momentum Settings)",
    os.path.join(outdir, "all_momenta_stability_pcal.png")
)

for mom in available_momenta:
    momdf = df[df["momentum"] == mom].copy()
    make_plot(
        momdf,
        f"SHMS Calorimeter Stability ({mom:g} GeV)",
        os.path.join(outdir, f"stability_pcal_{mom:g}GeV.png")
    )

print(f"Wrote {1 + len(available_momenta)} plots to {outdir}")
