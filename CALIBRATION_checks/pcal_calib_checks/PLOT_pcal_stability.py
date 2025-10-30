#!/usr/bin/env python3

import os, re
import pandas as pd
import matplotlib.pyplot as plt

fit_results_filepath = "DAT/FIT_pcal_results.dat"
auxfiles_runlist_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat"
processed_fit_filepath = "DAT/FIT_pcal_results_processed.dat"

runnums = []
shms_p = []
runtypes = []

def extract_parts(line):
    parts = re.split(r'\s+', line.strip())
    runnum = None
    momentum = None
    run_type = ""
    run_types = ["hole", "optics", "heep", "hee", "hmsdis", "shmsdis", "pi-sidis", "pi+sidis", "junk"]

    if len(parts) > 0:
        runnum = parts[0]
    if len(parts) >= 9:
        momentum_placeholder = parts[8]
        try:
            momentum = float(momentum_placeholder)
        except ValueError:
            m = re.match(r"^[+\-]?\d+(?:\.\d+)?", momentum_placeholder)
            if m:
                momentum = float(m.group(0))
    for part in parts:
        part_lower = part.strip().lower()
        if part_lower in run_types:
            run_type = part_lower
            break

    return runnum, momentum, run_type

        
with open(auxfiles_runlist_filepath, "r") as infile:
    lines = infile.readlines()

    for line in lines:
        if line.lstrip().startswith("#"):        # dropping comment lines
            continue
        if line.lstrip().startswith("!"):        # dropping lines starting with !
            continue
        if re.match(r"^\s*[-=*]", line):         # dropping separator lines
            continue

        runnum, momentum, run_type = extract_parts(line)
        if runnum is not None and momentum is not None and run_type is not None:
            runnums.append(runnum)
            shms_p.append(momentum)
            runtypes.append(run_type)

run_to_mom = dict(zip(runnums, shms_p))
run_to_type = dict(zip(runnums, runtypes))

with open(fit_results_filepath, "r") as infile, open(processed_fit_filepath, "w") as outfile:
    outfile.write("#runnum\tfit_mean\tmean_err\tfit_sigma\tsigma_err\tshms_p\truntype\n") # Writing the header
    for line in infile:
        if line.lstrip().startswith("#"): # And skipping here the header from the infile
            continue
        parts = re.split(r'\s+', line.strip())
        if len(parts) < 3:
            continue
        runnum = parts[0]
        fit_mean = parts[1]
        mean_err = parts[2]
        fit_sigma = parts[3]
        sigma_err = parts[4]
        momentum = run_to_mom.get(runnum, "")
        run_type = run_to_type.get(runnum, "")

        if momentum is None:
            raise ValueError(f"Run {runnum} has no momentum assigned!")
        if run_type is None or pd.isna(run_type):
            raise ValueError(f"Run {runnum} has no run type assigned!")
        outfile.write(f"{runnum}\t{fit_mean}\t{mean_err}\t{fit_sigma}\t{sigma_err}\t{momentum}\t{run_type}\n")

print(f"Wrote processed fit results to {processed_fit_filepath}")

df = pd.read_csv(
    processed_fit_filepath,
    comment="#",
    sep=r"\s+",
    names=["runnum", "fit_mean", "mean_err", "fit_sigma", "sigma_err", "momentum", "run_type"],
    dtype={"runnum": str, "momentum":str, "run_type":str}
)

df["fit_mean"] = df["fit_mean"].astype(float)
df["fit_sigma"] = df["fit_sigma"].astype(float)
df["runnum"] = df["runnum"].astype(str).str.replace(r"\.0$", "", regex=True)

df["runnum_int"] = df["runnum"].astype(int)
df = df.sort_values("runnum_int")

available_types = [str(rt) for rt in df["run_type"].unique() if isinstance(rt, str) or not pd.isna(rt)]
available_types = [rt for rt in available_types if rt != 'nan']
available_types = sorted(available_types)

print("Available run types with types:")
for i, val in enumerate(df["run_type"].unique()):
    print(i, repr(val), type(val))
print("\nAvailable run types:", ", ".join(available_types))

selected = input("Enter comma-separated run types to include (or 'all'): ").strip().lower()

if selected != "all":
    chosen_types = [s.strip() for s in selected.split(",")]
    df = df[df["run_type"].isin(chosen_types)]

if df.empty:
    print("No runs match the selected types.")
    exit()

df["runnum_int"] = df["runnum"].astype(int)
df = df.sort_values("runnum_int")

plt.figure(figsize=(10, 6))

for run_type, subdf in df.groupby("run_type"):
    plt.errorbar(
        subdf["runnum_int"].to_numpy(),
        subdf["fit_mean"].to_numpy(),
        yerr=subdf["mean_err"].to_numpy(),
        fmt="o",
        capsize=0,
        elinewidth = 0.5,
        markersize = 3,
        label=run_type
    )
plt.axhline(1.0, color='navy', linestyle='--', linewidth=1.2, label='y = 1')

plt.xlabel("Run Number", fontsize=13)
plt.ylabel("Fit Mean", fontsize=13)
plt.title(f"SHMS Calorimeter Fitted E/p vs Run Number\nSelected Run Types: {selected}", fontsize=15)
plt.legend(title="Run Type", fontsize=10)
plt.grid(True, linestyle="--", alpha=0.6)

if selected == "all":
    prefix = "all"
else:
    prefix = "_".join(chosen_types)
    prefix = prefix.replace(" ", "_").replace("/", "-")

outfile = f"STABILITY_plots/{prefix}_etottracknorm_vs_runnum_pcal.png"

plt.tight_layout()
plt.savefig(outfile, dpi = 300)
plt.close()
print(f"Saved plot to {outfile}")

#########################

plt.figure(figsize=(10, 6))

for run_type, subdf in df.groupby("run_type"):
    plt.errorbar(
        subdf["runnum_int"].to_numpy(),
        subdf["fit_sigma"].to_numpy(),
        yerr=subdf["sigma_err"].to_numpy(),
        fmt="o",
        capsize=0,
        elinewidth = 0.5,
        markersize = 3,
        label=run_type
    )

plt.xlabel("Run Number", fontsize=13)
plt.ylabel("Fit Mean", fontsize=13)
plt.title(f"SHMS Calorimeter Fitted E/p Sigma vs Run Number\nSelected Run Types: {selected}", fontsize=15)
plt.legend(title="Run Type", fontsize=10)
plt.grid(True, linestyle="--", alpha=0.6)

if selected == "all":
    prefix = "all"
else:
    prefix = "_".join(chosen_types)
    prefix = prefix.replace(" ", "_").replace("/", "-")

outfile = f"STABILITY_plots/{prefix}_etottracknorm_vs_runnum_sigma_pcal.png"

plt.tight_layout()
plt.savefig(outfile, dpi = 300)
plt.close()
print(f"Saved plot to {outfile}")
