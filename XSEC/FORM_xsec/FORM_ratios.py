#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import os, sys

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.jra_nprat import jra_nprat
from INIT.config import parse_run_type, parse_beam_pass, parse_target

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None
arg3 = sys.argv[3] if len(sys.argv) > 3 else None
arg4 = sys.argv[4] if len(sys.argv) > 4 else None

selected_run_type = parse_run_type(arg1)
selected_beam_pass, beam_prefix = parse_beam_pass(arg2)
num_abbrev, num_longname, num_shortname, num_A, num_Z = parse_target(arg3)
denom_abbrev, denom_longname, denom_shortname, denom_A, denom_Z = parse_target(arg4)

A_ratio = num_A / denom_A
num_N = num_A - num_Z
denom_N = denom_A - denom_Z

# -----------------------------------------------------
# Filepaths
# -----------------------------------------------------
num_csv_filepath = f"{num_abbrev.upper()}/XSEC_{selected_run_type}_{selected_beam_pass}pass_{num_abbrev}.csv"

denom_csv_filepath = f"{denom_abbrev.upper()}/XSEC_{selected_run_type}_{selected_beam_pass}pass_{denom_abbrev}.csv"

# -----------------------------------------------------
# Preparing Dataframes
# -----------------------------------------------------
df_num = pd.read_csv(num_csv_filepath)

df_denom = pd.read_csv(denom_csv_filepath)

merge_cols = ["eprime", "theta", "xbj", "q2", "w", "epsilon"]

num_rename = {col: f"{col}_num" for col in df_num.columns if col not in merge_cols}

denom_rename = {col: f"{col}_denom" for col in df_denom.columns if col not in merge_cols}

df_num = df_num.rename(columns=num_rename)

df_denom = df_denom.rename(columns=denom_rename)

df_merged = pd.merge(df_num, df_denom, on=merge_cols, how="inner")

f2rat = np.array([jra_nprat(x, q) for x, q in zip(df_merged["xbj"], df_merged["q2"])])

df_merged["f2rat"] = f2rat

df_merged["iso_corr"] = (0.5*(num_Z + num_N) * (1.0 + df_merged["f2rat"]) / (num_Z + num_N * df_merged["f2rat"]))

# Defining here the output ratios NOT normalized per nucleon

df_merged["raw_ratio"] = df_merged["xsec_exp_num"] / df_merged["xsec_exp_denom"]

df_merged["raw_ratio_err"] = df_merged["raw_ratio"] * np.sqrt((df_merged["xsec_exp_err_num"] / df_merged["xsec_exp_num"])**2 + (df_merged["xsec_exp_err_denom"] / df_merged["xsec_exp_denom"])**2)

df_merged["xsec_ratio"] = df_merged["raw_ratio"] * df_merged["iso_corr"]

df_merged["xsec_ratio_err"] = df_merged["raw_ratio_err"] * df_merged["iso_corr"]

# Normalizing per nucleon now,

df_merged["xsec_ratio_per_nucleon"] = df_merged["xsec_ratio"] / A_ratio

df_merged["xsec_ratio_per_nucleon_err"] = df_merged["xsec_ratio_err"] / A_ratio

# -----------------------------------------------------
# Save output csv
# -----------------------------------------------------
output_dir = "RATIOS"
os.makedirs(output_dir, exist_ok=True)

output_csv_filepath = f"{output_dir}/XSEC_RATIO_{selected_run_type}_{selected_beam_pass}pass_{num_abbrev}_to_{denom_abbrev}.csv"

df_merged.to_csv(output_csv_filepath, index=False)

print(f"Saved → {output_csv_filepath}")

# -----------------------------------------------------
# Plotting
# -----------------------------------------------------
output_pdf_filepath = f"{output_dir}/XSEC_RATIO_{selected_run_type}_{selected_beam_pass}pass_{num_abbrev}_to_{denom_abbrev}.pdf"

vars_to_plot = {
    "eprime": df_merged["eprime"].to_numpy(),
    "xbj": df_merged["xbj"].to_numpy(),
    "q2": df_merged["q2"].to_numpy(),
    "w": df_merged["w"].to_numpy(),
    "epsilon": df_merged["epsilon"].to_numpy(),
    }

xsec_ratio_per_nucleon = df_merged["xsec_ratio_per_nucleon"].to_numpy()
xsec_ratio_per_nucleon_err = df_merged["xsec_ratio_per_nucleon_err"].to_numpy()

pp = PdfPages(output_pdf_filepath)

for var, val in vars_to_plot.items():
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.errorbar(val, xsec_ratio_per_nucleon, yerr=xsec_ratio_per_nucleon_err, fmt='o', markersize=6, capsize=0, color='navy')
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.set_xlabel(f"{var}")
    ax.set_ylabel("Cross Section Ratio per Nucleon")
    ax.set_title(f"{selected_run_type.upper()} {selected_beam_pass}Pass {num_shortname}/{denom_shortname} Cross Section Ratio per Nucleon")
    ax.grid()

    pp.savefig(fig)
    plt.close(fig)
    
pp.close()
print(f"Saved PDF → {output_pdf_filepath}")


    
