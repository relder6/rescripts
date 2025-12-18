#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import os

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
selected_type = input("Enter desired run type (default HMSDIS): ").strip().lower()
if not selected_type:
    selected_type = "hmsdis"

selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ").strip()
selected_beam_pass_to_energy_prefix = {
    "1": "2.", "2": "4.", "3": "6.", "4": "8.", "5": "10."
}
beam_prefix = selected_beam_pass_to_energy_prefix.get(selected_beam_pass)
if not beam_prefix:
    print(f"Unknown pass: {selected_beam_pass}.  Please try again.")
    exit(1)

selected_target_shortcut_to_target_variable = {
    "al":"al","al13":"al","aluminum":"al",
    "c":"c","c12":"c","carbon":"c",
    "cu":"cu","cu29":"cu","copper":"cu",
    "opt1":"optics1","optics1":"optics1",
    "opt2":"optics2","optics2":"optics2",
    "d2":"ld2","ld2":"ld2",
    "h2":"lh2","lh2":"lh2",
    "hole":"hole","chole":"hole","c-hole":"hole",
    "dummy":"dummy","dum":"dummy"
}

selected_target_shortname_to_title_longname = {
    "al":"Aluminum",
    "c":"Carbon",
    "cu":"Copper",
    "opt1":"Optics1",
    "opt2":"Optics2",
    "ld2":"Deuterium",
    "lh2":"Hydrogen",
    "hole":"Carbon Hole",
    "dummy":"Dummy"}

# Targets
selected_num = input("Enter numerator target: ").strip().lower()
num_short = selected_target_shortcut_to_target_variable.get(selected_num)
if not num_short:
    print(f"Unknown target: {selected_num}. Please try again.")
    exit(1)
num_long = selected_target_shortname_to_title_longname[num_short]
    
selected_denom = input("Enter denominator target: ").strip().lower()
denom_short = selected_target_shortcut_to_target_variable.get(selected_denom)
if not denom_short:
    print(f"Unknown target: {selected_denom}. Please try again.")
    exit(1)
denom_long = selected_target_shortname_to_title_longname[denom_short]

# -----------------------------------------------------
# Filepaths
# -----------------------------------------------------
num_csv_filepath = f"{num_short.upper()}/XSEC_{selected_type}_{selected_beam_pass}pass_{num_short}.csv"

denom_csv_filepath = f"{denom_short.upper()}/XSEC_{selected_type}_{selected_beam_pass}pass_{denom_short}.csv"

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

# Ratio of experimental cross sections
df_merged["xsec_ratio"] = df_merged["xsec_exp_num"] / df_merged["xsec_exp_denom"]

# Error propagation for ratio: (deltaR / R)^2 = (deltaA / A)^2 + (deltaB / B)^2
df_merged["xsec_ratio_err"] = df_merged["xsec_ratio"] * np.sqrt((df_merged["xsec_exp_err_num"] / df_merged["xsec_exp_num"])**2 + (df_merged["xsec_exp_err_denom"] / df_merged["xsec_exp_denom"])**2)

# -----------------------------------------------------
# Save output csv
# -----------------------------------------------------
output_dir = "RATIOS"
os.makedirs(output_dir, exist_ok=True)

output_csv_filepath = f"{output_dir}/XSEC_RATIO_{selected_type}_{selected_beam_pass}pass_{num_short}_to_{denom_short}.csv"

df_merged.to_csv(output_csv_filepath, index=False)

print(f"Saved → {output_csv_filepath}")

# -----------------------------------------------------
# Plotting
# -----------------------------------------------------
output_pdf_filepath = f"{output_dir}/XSEC_RATIO_{selected_type}_{selected_beam_pass}pass_{num_short}_to_{denom_short}.pdf"

vars_to_plot = {
    "eprime": df_merged["eprime"].to_numpy(),
    "xbj": df_merged["xbj"].to_numpy(),
    "q2": df_merged["q2"].to_numpy(),
    "w": df_merged["w"].to_numpy(),
    "epsilon": df_merged["epsilon"].to_numpy(),
    }

xsec_ratio = df_merged["xsec_ratio"].to_numpy()
xsec_ratio_err = df_merged["xsec_ratio_err"].to_numpy()

pp = PdfPages(output_pdf_filepath)

for var, val in vars_to_plot.items():
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.errorbar(val, xsec_ratio, yerr=xsec_ratio_err, fmt='o', markersize=3, capsize=0)
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.set_xlabel(f"{var}")
    ax.set_ylabel("Cross Section Ratio")
    ax.set_title(f"{num_long} to {denom_long} {selected_type.upper()} Experimental Cross Section Ratio at {selected_beam_pass}Pass")
    ax.grid()

    pp.savefig(fig)
    plt.close(fig)
    
pp.close()
print(f"Saved PDF → {output_pdf_filepath}")


    
