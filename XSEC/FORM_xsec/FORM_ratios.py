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


# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
if len(sys.argv) == 5:
    selected_run_type = sys.argv[1].strip().lower()
    selected_beam_pass = sys.argv[2].strip()
    selected_num = sys.argv[3].strip().lower()
    selected_denom = sys.argv[4].strip().lower()
else:
    selected_run_type = input("Enter desired run type (default HMSDIS): ").strip().lower()
    if not selected_run_type:
        selected_run_type = "hmsdis"
    selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ").strip()
    selected_num = input("Enter numerator target: ").strip().lower()
    selected_denom = input("Enter denominator target: ").strip().lower()
    
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
num_short = selected_target_shortcut_to_target_variable.get(selected_num)
if not num_short:
    print(f"Unknown target: {selected_num}. Please try again.")
    exit(1)
num_long = selected_target_shortname_to_title_longname[num_short]
    
denom_short = selected_target_shortcut_to_target_variable.get(selected_denom)
if not denom_short:
    print(f"Unknown target: {selected_denom}. Please try again.")
    exit(1)
denom_long = selected_target_shortname_to_title_longname[denom_short]

selected_target_shortname_to_AZ = {"al":   (27, 13),
                                   "c":    (12, 6),
                                   "cu":   (64, 29),
                                   "ld2":  (2, 1),
                                   "lh2":  (1, 1)}

A_num, Z_num = selected_target_shortname_to_AZ.get(num_short)

N_num = A_num - Z_num

A_denom, Z_denom = selected_target_shortname_to_AZ.get(denom_short, (0.0, 0.0))

A_ratio = A_num / A_denom

# -----------------------------------------------------
# Filepaths
# -----------------------------------------------------
num_csv_filepath = f"{num_short.upper()}/XSEC_{selected_run_type}_{selected_beam_pass}pass_{num_short}.csv"

denom_csv_filepath = f"{denom_short.upper()}/XSEC_{selected_run_type}_{selected_beam_pass}pass_{denom_short}.csv"

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

df_merged["xsec_ratio"] = df_merged["xsec_exp_num"] / df_merged["xsec_exp_denom"] 

df_merged["xsec_ratio_err"] = df_merged["xsec_ratio"] * np.sqrt((df_merged["xsec_exp_err_num"] / df_merged["xsec_exp_num"])**2 + (df_merged["xsec_exp_err_denom"] / df_merged["xsec_exp_denom"])**2)

# Normalizing per nucleon

df_merged["xsec_ratio_norm"] = df_merged["xsec_ratio"] / A_ratio

df_merged["xsec_ratio_norm_err"] = df_merged["xsec_ratio_err"] / A_ratio

# Isoscalar corrections

f2rat = np.array([jra_nprat(x, q) for x, q in zip(df_merged["xbj"], df_merged["q2"])])

df_merged["f2rat"] = f2rat

df_merged["iso_corr"] = (0.5*(Z_num + N_num) * (1.0 + df_merged["f2rat"]) / (Z_num + N_num * df_merged["f2rat"]))

df_merged["xsec_ratio_final"] = df_merged["xsec_ratio_norm"] * df_merged["iso_corr"]

df_merged["xsec_ratio_final_err"] = df_merged["xsec_ratio_norm_err"] * df_merged["iso_corr"]

# -----------------------------------------------------
# Save output csv
# -----------------------------------------------------
output_dir = "RATIOS"
os.makedirs(output_dir, exist_ok=True)

output_csv_filepath = f"{output_dir}/XSEC_RATIO_{selected_run_type}_{selected_beam_pass}pass_{num_short}_to_{denom_short}.csv"

df_merged.to_csv(output_csv_filepath, index=False)

print(f"Saved → {output_csv_filepath}")

# -----------------------------------------------------
# Plotting
# -----------------------------------------------------
output_pdf_filepath = f"{output_dir}/XSEC_RATIO_{selected_run_type}_{selected_beam_pass}pass_{num_short}_to_{denom_short}.pdf"

vars_to_plot = {
    "eprime": df_merged["eprime"].to_numpy(),
    "xbj": df_merged["xbj"].to_numpy(),
    "q2": df_merged["q2"].to_numpy(),
    "w": df_merged["w"].to_numpy(),
    "epsilon": df_merged["epsilon"].to_numpy(),
    }

xsec_ratio_norm = df_merged["xsec_ratio_norm"].to_numpy()
xsec_ratio_norm_err = df_merged["xsec_ratio_norm_err"].to_numpy()

pp = PdfPages(output_pdf_filepath)

for var, val in vars_to_plot.items():
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.errorbar(val, xsec_ratio_norm, yerr=xsec_ratio_norm_err, fmt='o', markersize=6, capsize=0, color='navy')
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.set_xlabel(f"{var}")
    ax.set_ylabel("Cross Section Ratio per Nucleon")
    ax.set_title(f"{selected_run_type.upper()} {selected_beam_pass}Pass {num_long}/{denom_long} Cross Section Ratio per Nucleon")
    ax.grid()

    pp.savefig(fig)
    plt.close(fig)
    
pp.close()
print(f"Saved PDF → {output_pdf_filepath}")


    
