#!/usr/bin/env python3

import os, re, sys, csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fit_results_filepath = "DAT/FIT_hcal_results.dat"

rows = []
with open(fit_results_filepath, "r") as f:
    reader = csv.DictReader(filter(lambda r: r.strip() and not r.startswith("#"), f), delimiter="\t")
    for r in reader:
        try:
            rows.append({
                "runnum": str(r["runnum"]).strip(),
                "run_type": str(r["runtype"]).strip().lower(),
                "fit_mean": float(r["fit_mean"]),
                "mean_err": float(r["mean_err"]),
                "fit_sigma": float(r["fit_sigma"]),
                "sigma_err": float(r["sigma_err"])
            })
        except (ValueError, TypeError):
            continue

df = pd.DataFrame(rows)
if df.empty:
    print("No valid fit rows found.")
    sys.exit(1)

df["runnum"] = df["runnum"].str.replace(r"\.0$", "", regex=True)
df["runnum_int"] = df["runnum"].astype(int)
df = df.sort_values("runnum_int")

available_types = sorted(df["run_type"].unique())

if len(sys.argv) > 1:
    selected = sys.argv[1].strip().lower()
else:
    print("\nAvailable run types:", ", ".join(available_types))
    selected = input("Enter comma-separated run types to include (or 'all'): ").strip().lower()

if selected != "all":
    chosen_types = [s.strip() for s in selected.split(",")]
    df = df[df["run_type"].isin(chosen_types)]

if df.empty:
    print("No runs match the selected types.")
    sys.exit(0)

df["residual"] = df["fit_mean"] - 1.0

fig, (ax_top, ax_bot) = plt.subplots(nrows=2, ncols=1, figsize=(8,8), gridspec_kw={"height_ratios":[3,1]})

for run_type, subdf in df.groupby("run_type"):
    ax_top.errorbar(subdf["runnum_int"].to_numpy(), subdf["fit_mean"].to_numpy(),
                    yerr=subdf["mean_err"].to_numpy(),
                    fmt="o", capsize=0, elinewidth=0.5, markersize=3, label=run_type)

ax_top.axhline(1.0, color="navy", linestyle="--", linewidth=1.2, label="y = 1")
ax_top.set_ylabel(r"$\mu_{\rm fit}$, H.cal.etottracknorm", fontsize=13)
ax_top.set_xlabel("Run Number", fontsize=13)
ax_top.legend(title="Run Type", fontsize=10)
ax_top.grid(True, linestyle="--", alpha=0.6)

ax_bot.hist(df["residual"], bins=100, histtype="stepfilled", alpha=0.7, edgecolor="black")
ax_bot.axvline(0, color="black", linestyle="--", linewidth=1)

sigma = df["residual"].std()
for n in [1,2,3]:
    ax_bot.axvline(n*sigma, color="red", linestyle="--", alpha=0.6)
    ax_bot.axvline(-n*sigma, color="red", linestyle="--", alpha=0.6)

ax_bot.set_xlabel("Residual (Fit Mean − 1.0)", fontsize=13)
ax_bot.set_ylabel("Count", fontsize=13)
ax_bot.grid(True, linestyle="--", alpha=0.6)

prefix = "all" if selected == "all" else "_".join(chosen_types).replace(" ", "_").replace("/", "-")
outfile_mean = f"STABILITY_plots/{prefix}_stability_hcal.png"

fig.suptitle(f"HMS Calorimeter Fitted E/p Mean vs Run Number and Residuals\nSelected Run Types: {selected}", fontsize=15)
plt.tight_layout()
plt.savefig(outfile_mean, dpi=300)
plt.close()

print(f"Saved plot to {outfile_mean}")
