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

selected_target = input("Enter desired target (options: C, Cu, Al, LD2, LH2, Dummy): ").strip().lower()
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

selected_target_shortname = selected_target_shortcut_to_target_variable.get(selected_target)
if not selected_target_shortname:
    print(f"Unknown target: {selected_target}.  Please try again.")
    exit(1)
    
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

selected_target_title_longname = selected_target_shortname_to_title_longname.get(selected_target_shortname)

# -----------------------------------------------------
# Filepaths
# -----------------------------------------------------
model_xsec_filepath = f"../MODEL_xsec/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_model_xsec.csv"

data_to_mc_filepath = f"../DATA_to_MC/{selected_target_shortname.upper()}/DATA_to_MC_{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_H_gtr_p.csv"

xsec_pdf_output = f"PDFs/XSEC_{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}.pdf"

# -----------------------------------------------------
# Preparing Dataframes
# -----------------------------------------------------
df_model_xsec = pd.read_csv(model_xsec_filepath)

df_data_to_mc = pd.read_csv(data_to_mc_filepath).rename(columns={"bin_center": "eprime"})

# -----------------------------------------------------
# Forming the xsec
# -----------------------------------------------------
df_model_xsec["modelxsec"] = pd.to_numeric(df_model_xsec["modelxsec"], errors="coerce")

df_merged = pd.merge(df_data_to_mc, df_model_xsec, on="eprime", how="inner")

df_merged["xsec_exp"] = df_merged["ratio"] * df_merged["modelxsec"]

# -----------------------------------------------------
# Error propagation
# -----------------------------------------------------
df_merged["xsec_exp_err"] = df_merged["xsec_exp"] * (df_merged["ratio_err"] / df_merged["ratio"])

df_merged["xsec_exp"] = df_merged["xsec_exp"].replace([np.inf, -np.inf], np.nan)
df_merged["xsec_exp_err"] = df_merged["xsec_exp_err"].replace([np.inf, -np.inf], np.nan)

# -----------------------------------------------------
# Save output csv
# -----------------------------------------------------
output_dir = f"{selected_target_shortname.upper()}"
os.makedirs(output_dir, exist_ok=True)

output_filepath = f"{output_dir}/XSEC_{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}.csv"

final_columns = ["eprime", "theta", "xbj", "q2", "w", "modelxsec", "xsec_exp", "xsec_exp_err"]

df_final = df_merged[final_columns]

df_final.to_csv(output_filepath, index=False)

print(f"Saved → {output_filepath}")

# -----------------------------------------------------
# Plotting
# -----------------------------------------------------
vars_to_plot = {
    "eprime": df_final["eprime"].to_numpy(),
    "xbj": df_final["xbj"].to_numpy(),
    "q2": df_final["q2"].to_numpy(),
    "w": df_final["w"].to_numpy(),
}

xsec_final = df_final["xsec_exp"].to_numpy()
xsec_err_final = df_final["xsec_exp_err"].to_numpy()

pp = PdfPages(xsec_pdf_output)

for var, val in vars_to_plot.items():
    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    mask = ~np.isnan(xsec_final)
    model_masked = df_final["modelxsec"].to_numpy().copy()
    model_masked[~mask] = np.nan

    
    ax.errorbar(val, xsec_final, yerr=xsec_err_final, fmt='o', markersize=3, capsize=0, label = "Data")
    ax.plot(val, model_masked, linestyle = ":", marker = "", label = "Model")
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.set_xlabel(f"{var}")
    ax.set_ylabel("Cross Section (ub/GeV/sr)")
    ax.set_title(f"{selected_target_title_longname} {selected_type.upper()} Experimental Cross Section at {selected_beam_pass}Pass")
    ax.grid()
    ax.legend()
    ymax = np.nanmax([xsec_final, df_final["modelxsec"].to_numpy()])
    ax.set_ylim(0, ymax*1.1)
    pp.savefig(fig)
    plt.close(fig)

pp.close()
print(f"Saved PDF → {xsec_pdf_output}")
