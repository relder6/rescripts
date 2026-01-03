#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import os, re, sys
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)    
from INIT import get_common_run_inputs

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
# selected_type = input("Enter desired run type (default HMSDIS): ").strip().lower()
# if not selected_type:
#     selected_type = "hmsdis"

# selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ").strip()
# selected_beam_pass_to_energy_prefix = {
#     "1": "2.", "2": "4.", "3": "6.", "4": "8.", "5": "10."
# }
# beam_prefix = selected_beam_pass_to_energy_prefix.get(selected_beam_pass)
# if not beam_prefix:
#     print(f"Unknown pass: {selected_beam_pass}.  Please try again.")
#     exit(1)

# selected_target = input("Enter desired target (options: C, Cu, Al, LD2, LH2, Dummy): ").strip().lower()
# selected_target_shortcut_to_target_variable = {
#     "al":"al","al13":"al","aluminum":"al",
#     "c":"c","c12":"c","carbon":"c",
#     "cu":"cu","cu29":"cu","copper":"cu",
#     "opt1":"optics1","optics1":"optics1",
#     "opt2":"optics2","optics2":"optics2",
#     "d2":"ld2","ld2":"ld2",
#     "h2":"lh2","lh2":"lh2",
#     "hole":"hole","chole":"hole","c-hole":"hole",
#     "dummy":"dummy","dum":"dummy"
# }

# selected_target_shortname = selected_target_shortcut_to_target_variable.get(selected_target)
# if not selected_target_shortname:
#     print(f"Unknown target: {selected_target}.  Please try again.")
#     exit(1)
    
# selected_target_shortname_to_title_longname = {
#     "al":"Aluminum",
#     "c":"Carbon",
#     "cu":"Copper",
#     "opt1":"Optics1",
#     "opt2":"Optics2",
#     "ld2":"Deuterium",
#     "lh2":"Hydrogen",
#     "hole":"Carbon Hole",
#     "dummy":"Dummy"}

# selected_target_title_longname = selected_target_shortname_to_title_longname.get(selected_target_shortname)
# =====================================================================
# Handling user inputs
# =====================================================================
selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, selected_target_titlename = get_common_run_inputs()


# -----------------------------------------------------
# Files
# -----------------------------------------------------
base_dir = f"../MAKE_csvs/{selected_target_shortname.upper()}"
    
variables = ["H_gtr_dp", "H_gtr_ph", "H_gtr_th", "H_kin_Q2", "H_kin_x_bj", "H_kin_W", "H_gtr_p"]

pdf_output = f"PDFs/DATA_to_MC_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}.pdf"
pp = PdfPages(pdf_output)

# -----------------------------------------------------
# Loop over variables and produce figures for SOLID targets
# -----------------------------------------------------
if selected_target_shortname in {"al", "c", "cu"}:
    for var in variables:
        histo_filepath = f"{base_dir}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_histo.csv"
        err_filepath   = f"{base_dir}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_err.csv"

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
        # Saving binned yield ratios to output CSV
        # --------------------------------------------------
        output_dir = f"{selected_target_shortname.upper()}"
        output_filepath = f"{output_dir}/DATA_to_MC_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}.csv"
        os.makedirs(output_dir, exist_ok=True)

        df_output = pd.DataFrame({
            "bin_center": np.round(bin_centers, 4),
            "Y_sub": Y_sub,
            "E_sub": E_sub,
            "Y_mc": Y_mc,
            "E_mc": E_mc,
            "ratio": ratio,
            "ratio_err": ratio_err
        })

        df_output.to_csv(output_filepath, index=False)
        print(f"Saved binned yield CSV → {output_filepath}")
        
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
        ax_top.set_title(f"{selected_target_titlename} {var}: Data (e⁻ - e⁺) to MC Yields")
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
    exit(1)

# -----------------------------------------------------
# Loop over variables and produce figures for LIQUID targets
# -----------------------------------------------------
if selected_target_shortname in {"ld2", "lh2"}:
    dummy_dir = f"../MAKE_csvs/DUMMY"

    R_dummy = 1 / 3.7825 if selected_target_shortname == 'ld2' else 1 / 3.5274

    # Ensure H_gtr_dp is first
    ordered_variables = ["H_gtr_dp"] + [v for v in variables if v != "H_gtr_dp"]

    for var in ordered_variables:
        # -------------------------
        # Read CSVs
        # -------------------------
        histo_filepath = f"{base_dir}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_histo.csv"
        err_filepath   = f"{base_dir}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_err.csv"

        df_yield = pd.read_csv(histo_filepath)
        df_err   = pd.read_csv(err_filepath)

        runinfo_cols = ['runnum','charge','polarity']
        bin_cols = [c for c in df_yield.columns if c not in runinfo_cols]

        # Dummy
        dummy_histo_filepath = f"{dummy_dir}/{selected_run_type}_{selected_beam_pass}pass_dummy_{var}_histo.csv"
        dummy_err_filepath   = f"{dummy_dir}/{selected_run_type}_{selected_beam_pass}pass_dummy_{var}_err.csv"

        df_yield_dummy = pd.read_csv(dummy_histo_filepath)
        df_err_dummy   = pd.read_csv(dummy_err_filepath)
        bin_cols_dummy = [c for c in df_yield_dummy.columns if c not in runinfo_cols]

        # -------------------------
        # Masks
        # -------------------------
        elec_mask = df_yield['polarity'] == '-'
        pos_mask  = df_yield['polarity'] == '+'
        mc_mask   = df_yield['polarity'] == 'mc'

        elec_mask_dummy = df_yield_dummy['polarity'] == '-'
        pos_mask_dummy  = df_yield_dummy['polarity'] == '+'

        # -------------------------
        # Charge normalization
        # -------------------------
        charge_norm_elec = 1000 / df_yield.loc[elec_mask, 'charge'].sum()
        charge_norm_pos  = 1000 / df_yield.loc[pos_mask , 'charge'].sum()
        df_yield.loc[elec_mask, bin_cols] *= charge_norm_elec
        df_yield.loc[pos_mask , bin_cols] *= charge_norm_pos
        df_err.loc[elec_mask, bin_cols]   *= charge_norm_elec
        df_err.loc[pos_mask , bin_cols]   *= charge_norm_pos

        charge_norm_elec_dummy = 1000 / df_yield_dummy.loc[elec_mask_dummy, 'charge'].sum()
        charge_norm_pos_dummy  = 1000 / df_yield_dummy.loc[pos_mask_dummy, 'charge'].sum()
        df_yield_dummy.loc[elec_mask_dummy, bin_cols_dummy] *= charge_norm_elec_dummy
        df_yield_dummy.loc[pos_mask_dummy , bin_cols_dummy] *= charge_norm_pos_dummy
        df_err_dummy.loc[elec_mask_dummy, bin_cols_dummy] *= charge_norm_elec_dummy
        df_err_dummy.loc[pos_mask_dummy , bin_cols_dummy] *= charge_norm_pos_dummy

        # -------------------------
        # Sum yields
        # -------------------------
        Y_elec = df_yield.loc[elec_mask, bin_cols].sum(axis=0).values
        Y_pos  = df_yield.loc[pos_mask , bin_cols].sum(axis=0).values
        Y_mc   = df_yield.loc[mc_mask , bin_cols].values[0]

        Y_elec_dummy = df_yield_dummy.loc[elec_mask_dummy, bin_cols_dummy].sum(axis=0).values
        Y_pos_dummy  = df_yield_dummy.loc[pos_mask_dummy , bin_cols_dummy].sum(axis=0).values

        # -------------------------
        # Errors
        # -------------------------
        E_elec = np.sqrt(np.sum(df_err.loc[elec_mask, bin_cols].values**2, axis=0))
        E_pos  = np.sqrt(np.sum(df_err.loc[pos_mask , bin_cols].values**2, axis=0))
        E_mc   = df_err.loc[mc_mask, bin_cols].values[0]

        E_elec_dummy = np.sqrt(np.sum(df_err_dummy.loc[elec_mask_dummy, bin_cols_dummy].values**2, axis=0))
        E_pos_dummy  = np.sqrt(np.sum(df_err_dummy.loc[pos_mask_dummy , bin_cols_dummy].values**2, axis=0))

        # -------------------------
        # Dummy subtraction
        # -------------------------
        Y_sub = Y_elec - Y_pos - R_dummy * (Y_elec_dummy - Y_pos_dummy)
        E_sub = np.sqrt(E_elec**2 + E_pos**2 + (R_dummy * E_elec_dummy)**2 + (R_dummy * E_pos_dummy)**2)

        # Totals and positron fraction
        Y_elec_tot = np.nansum(Y_elec)
        Y_pos_tot  = np.nansum(Y_pos)
        Y_sub_tot  = np.nansum(Y_sub)
        pos_frac   = 100.0 * Y_pos_tot / Y_elec_tot if Y_elec_tot > 0 else np.nan

        # -------------------------
        # Ratio
        # -------------------------
        ratio = np.divide(Y_sub, Y_mc, out=np.full_like(Y_sub, np.nan), where=(Y_mc>0))
        ratio_err = np.full_like(ratio, np.nan)
        good = (Y_sub>0) & (Y_mc>0)
        ratio_err[good] = ratio[good] * np.sqrt((E_sub[good]/Y_sub[good])**2 + (E_mc[good]/Y_mc[good])**2)

        # -------------------------
        # Bin centers
        # -------------------------
        bin_left_edges = np.array([float(c) for c in bin_cols])
        bin_width = bin_left_edges[1] - bin_left_edges[0]
        bin_centers = bin_left_edges + 0.5 * bin_width

        # -------------------------
        # Save CSV
        # -------------------------
        output_dir = f"{selected_target_shortname.upper()}"
        os.makedirs(output_dir, exist_ok=True)
        df_output = pd.DataFrame({
            "bin_center": np.round(bin_centers,4),
            "Y_sub": Y_sub,
            "E_sub": E_sub,
            "Y_mc": Y_mc,
            "E_mc": E_mc,
            "ratio": ratio,
            "ratio_err": ratio_err
        })
        df_output.to_csv(f"{output_dir}/DATA_to_MC_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}.csv", index=False)
        print(f"Saved binned yield CSV → {output_dir}/DATA_to_MC_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}.csv")

        # -------------------------
        # Plotting
        # -------------------------
        fig, (ax_top, ax_bot) = plt.subplots(2,1,figsize=(8.5,11/2), gridspec_kw={'height_ratios':[3,1]}, sharex=True)

        # Always plot main data and MC scatter points
        ax_top.errorbar(bin_centers, Y_sub, yerr=E_sub, fmt='o', markersize=3, color='navy', label="Data")
        ax_top.errorbar(bin_centers, Y_mc,  yerr=E_mc, fmt='o', markersize=3, color='red',  label="MC")

        # Enhanced info for H_gtr_dp only
        if var == "H_gtr_dp":
            ax_top.errorbar(bin_centers, Y_elec, yerr=E_elec, fmt='o', markersize=3, color='darkorange', label='e⁻ yield')
            ax_top.errorbar(bin_centers, Y_pos,  yerr=E_pos,  fmt='o', markersize=3, color='darkgreen',  label='e⁺ yield')
            
            textstr = (
                f"Total e⁻ yield: {Y_elec_tot:8.2f}\n"
                f"Total e⁺ yield: {Y_pos_tot:8.2f}\n"
                f"Data (e⁻ - e⁺) :     {Y_sub_tot:8.2f}\n"
                f"Total MC:       {np.nansum(Y_mc):8.2f}\n"
                f"e⁺ fraction:    {pos_frac:6.2f} %\n"
                f"Data/MC ratio:  {100.0 * np.nansum(Y_sub)/np.nansum(Y_mc):6.2f} %"
            )

            ax_top.text(0.02, 0.75, textstr, transform=ax_top.transAxes, fontsize=10,
                        verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", alpha=0.9))
            
            plt.tight_layout()

        # Bottom panel: ratio
        ax_bot.errorbar(bin_centers, ratio, yerr=ratio_err, fmt='o', markersize=3, color='darkmagenta')
        ax_bot.axhline(1, color='gray', linestyle='--')
        ax_bot.set_ylabel("Data/MC"); ax_bot.set_ylim(0.75,1.25)
        ax_bot.set_xlabel(var); ax_bot.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        ax_bot.grid()

        ax_top.set_ylabel("Charge-Normalized Weighted Counts")
        ax_top.set_title(f"{selected_target_titlename} {var}: Data (e⁻ - e⁺) to MC Yields")
        ax_top.grid(); ax_top.legend(loc='upper right')

        fig.tight_layout(); fig.subplots_adjust(top=0.93, bottom=0.12, hspace=0.05)
        pp.savefig(fig); plt.close(fig)

    pp.close()
    print(f"\nSaved PDF → {pdf_output}")
    exit(1)


# -----------------------------------------------------
# Exception for other types (HEE, HEEP, Optics, HOLE)
# -----------------------------------------------------
else:
    print(f"\nTarget {selected_target_shortname} not currently supported by this script.")
    exit(1)
