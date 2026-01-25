#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
from matplotlib.ticker import PercentFormatter
import os, re, sys
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import get_common_run_inputs

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, selected_target_titlename, selected_target_A, selected_target_Z = get_common_run_inputs()

# -----------------------------------------------------
# Filepaths
# -----------------------------------------------------
model_xsec_filepath = f"../MODEL_xsec/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_model_xsec.csv"

data_to_mc_filepath = f"../DATA_to_MC/{selected_target_shortname.upper()}/DATA_to_MC_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_H_gtr_p.csv"

xsec_pdf_output = f"PDFs/XSEC_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}.pdf"

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

df_merged["A"] = selected_target_A
df_merged["Z"] = selected_target_Z

output_dir = f"{selected_target_shortname.upper()}"
os.makedirs(output_dir, exist_ok=True)

output_filepath = f"{output_dir}/XSEC_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}.csv"

final_columns = ["A", "Z", "eprime", "theta", "xbj", "q2", "w", "epsilon", "modelxsec", "xsec_exp", "xsec_exp_err"]

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
    "epsilon": df_final["epsilon"].to_numpy()
}

xsec_final = df_final["xsec_exp"].to_numpy()
xsec_err_final = df_final["xsec_exp_err"].to_numpy()

pp = PdfPages(xsec_pdf_output)

for var, val in vars_to_plot.items():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 7.5),
        gridspec_kw={"height_ratios":[3,1], "hspace":0.05}, sharex=True)

    mask = ~np.isnan(xsec_final)
    model = df_final["modelxsec"].to_numpy()

    ax1.errorbar(val[mask], xsec_final[mask], yerr=xsec_err_final[mask], fmt='o', markersize=3, color='navy', capsize=0, label="Data")
    ax1.plot(val[mask], model[mask], linestyle=":", color='red', label="Model")

    ax1.set_ylabel("Cross Section (μb/GeV/sr)")
    ax1.set_title(f"{selected_run_type.upper()} {selected_beam_pass}Pass {selected_target_titlename} Nuclear Cross Section")
    ax1.grid()
    ax1.legend()

    ymax = np.nanmax([xsec_final, model])
    ax1.set_ylim(0, ymax*1.1)

    residual = (xsec_final - model) / model
    residual_err = xsec_err_final / model

    ax2.errorbar(val[mask], residual[mask], yerr=residual_err[mask], fmt='o', markersize=3, capsize=0, color ='darkmagenta')
    ax2.axhline(0.0, color='k', linestyle='--', linewidth=1)
    ax2.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    ax2.set_xlabel(var)
    ax2.set_ylabel("Res. %")
    ax2.grid()

    ax2.relim()
    ax2.autoscale_view()

    pp.savefig(fig)
    plt.close(fig)

pp.close()

print(f"Saved PDF → {xsec_pdf_output}")
