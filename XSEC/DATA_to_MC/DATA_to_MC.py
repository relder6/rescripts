#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker

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
# Files
# -----------------------------------------------------
base_dir = f"../MAKE_csvs/{selected_target_shortname.upper()}"
if selected_target_shortname in {"ld2", "lh2"}:
    dummy_dir = f"../MAKE_csvs/DUMMY"
variables = ["H_gtr_dp", "H_gtr_ph", "H_gtr_th", "H_kin_Q2", "H_kin_x_bj", "H_kin_W"]

pdf_output = f"PDFs/DATA_to_MC_{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}.pdf"
pp = PdfPages(pdf_output)

# -----------------------------------------------------
# Loop over variables and produce figures
# -----------------------------------------------------
for var in variables:
    histo_filepath = f"{base_dir}/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_histo.csv"
    err_filepath   = f"{base_dir}/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_err.csv"

    df_yield = pd.read_csv(histo_filepath)
    df_err   = pd.read_csv(err_filepath)

    runinfo_cols = ['runnum','charge','polarity']
    bin_cols = [c for c in df_yield.columns if c not in runinfo_cols]

    # masks
    elec_mask = df_yield['polarity'] == '-'
    pos_mask  = df_yield['polarity'] == '+'
    mc_mask   = df_yield['polarity'] == 'mc'
    # --------------------------------------------------
    # CHARGE NORMALIZATION: scale data only so sum(charge)=1000
    # --------------------------------------------------
    total_charge_elec = df_yield.loc[elec_mask, 'charge'].sum()
    total_charge_pos = df_yield.loc[pos_mask, 'charge'].sum()
    
    charge_norm_elec = 1000 / total_charge_elec
    charge_norm_pos = 1000 / total_charge_pos

    df_yield.loc[elec_mask, bin_cols] *= charge_norm_elec
    df_err.loc[elec_mask, bin_cols]   *= charge_norm_elec

    df_yield.loc[pos_mask, bin_cols] *= charge_norm_pos
    df_err.loc[pos_mask, bin_cols] *= charge_norm_pos

    # --------------------------------------------------
    # Sum Yields (already normalized)
    # --------------------------------------------------
    Y_elec = df_yield.loc[elec_mask, bin_cols].astype(float).sum(axis=0).values
    Y_pos  = df_yield.loc[pos_mask , bin_cols].astype(float).sum(axis=0).values
    Y_mc   = df_yield.loc[mc_mask ,  bin_cols].astype(float).values[0]

    # --------------------------------------------------
    # Combine Errors
    # --------------------------------------------------
    E_elec = np.sqrt(np.sum(df_err.loc[elec_mask, bin_cols].astype(float).values**2, axis=0))
    E_pos  = np.sqrt(np.sum(df_err.loc[pos_mask , bin_cols].astype(float).values**2, axis=0))
    E_mc   = df_err.loc[mc_mask, bin_cols].astype(float).values[0]

    # --------------------------------------------------
    # Subtracted yield (e− − e+)
    # --------------------------------------------------
    Y_sub = Y_elec - Y_pos
    E_sub = np.sqrt(E_elec**2 + E_pos**2)

    # --------------------------------------------------
    # Ratio = (e− − e+) / MC
    # --------------------------------------------------
    ratio = np.divide(Y_sub, Y_mc, out=np.full_like(Y_sub, np.nan), where=(Y_mc > 0))
    ratio_err = np.full_like(ratio, np.nan)

    good = (Y_sub > 0) & (Y_mc > 0)
    ratio_err[good] = ratio[good] * np.sqrt(
        (E_sub[good] / Y_sub[good])**2 + (E_mc[good] / Y_mc[good])**2
    )

    # --------------------------------------------------
    # Bin centers
    # --------------------------------------------------
    bin_left_edges = np.array([float(cols) for cols in bin_cols])
    bin_width = bin_left_edges[1] - bin_left_edges[0]
    bin_centers = bin_left_edges + 0.5 * bin_width

    # --------------------------------------------------
    # Plotting
    # --------------------------------------------------
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.5, 11/2),
        gridspec_kw={'height_ratios':[3,1]},
        sharex=True
    )

    # --- top panel ---
    ax_top.errorbar(bin_centers, Y_sub, yerr=E_sub, fmt='o', markersize=3, capsize=0, color='navy', label="Data")
    ax_top.errorbar(bin_centers, Y_mc, yerr=E_mc, fmt='o', markersize=3, capsize=0, color='red', label="MC")
    ax_top.set_ylabel("Charge-Normalized Weighted Counts")
    ax_top.set_title(f"{selected_target_title_longname} {var}: Data (e⁻ - e⁺) to MC Yields")
    ax_top.grid()
    ax_top.legend()

    # --- bottom panel ---
    ax_bot.errorbar(bin_centers, ratio, yerr=ratio_err, fmt='o', markersize=3, capsize=0, color='darkmagenta')
    ax_bot.axhline(1, color='gray', linestyle='--')
    ax_bot.set_ylabel("Data/MC")
    ax_bot.set_ylim(0.75, 1.25)
    ax_bot.set_xlabel(var)
    ax_bot.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    ax_bot.grid()

    fig.tight_layout()
    fig.subplots_adjust(top = 0.93, bottom = 0.12, hspace = 0.05)
    
    pp.savefig(fig)
    plt.close(fig)

pp.close()
print(f"\nSaved PDF → {pdf_output}")
