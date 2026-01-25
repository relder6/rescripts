#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
import os, re, sys
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import get_common_run_inputs

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
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

        # Masks
        elec_mask = df_yield['polarity'] == '-'
        pos_mask  = df_yield['polarity'] == '+'
        mc_mask   = df_yield['polarity'] == 'mc'
        
        # Charge Normalization
        charge_norm_elec = 0.0
        if np.nansum(df_yield.loc[elec_mask, 'charge']) > 0:
            charge_norm_elec = 1000 / df_yield.loc[elec_mask, 'charge'].sum()
            df_yield.loc[elec_mask, bin_cols] *= charge_norm_elec
            df_err.loc[elec_mask, bin_cols] *= charge_norm_elec
        
        charge_norm_pos = 0.0
        if np.nansum(df_yield.loc[pos_mask, 'charge']) > 0:
            charge_norm_pos  = 1000 / df_yield.loc[pos_mask , 'charge'].sum()
            df_yield.loc[pos_mask , bin_cols] *= charge_norm_pos
            df_err.loc[pos_mask , bin_cols] *= charge_norm_pos
        
        # Sum Normalized Yields
        Yield_elec = np.zeros(len(bin_cols))
        if np.nansum(df_yield.loc[elec_mask, bin_cols]) > 0:
            Yield_elec = df_yield.loc[elec_mask, bin_cols].astype(float).sum(axis=0).values
            
        Yield_pos = np.zeros(len(bin_cols))
        if np.nansum(df_yield.loc[pos_mask, bin_cols]) > 0:
            Yield_pos  = df_yield.loc[pos_mask , bin_cols].astype(float).sum(axis=0).values

        Yield_mc = np.zeros(len(bin_cols))
        if np.nansum(df_yield.loc[mc_mask, bin_cols]) > 0:
            Yield_mc   = df_yield.loc[mc_mask ,  bin_cols].astype(float).values.flatten()

        # Error Propagation
        Err_elec = np.zeros(len(bin_cols))
        if np.nansum(df_err.loc[elec_mask, bin_cols]) > 0:
            Err_elec = np.sqrt(np.sum(df_err.loc[elec_mask, bin_cols].astype(float).values**2, axis=0))
            
        Err_pos = np.zeros(len(bin_cols))
        if np.nansum(df_err.loc[pos_mask, bin_cols]) > 0:    
            Err_pos  = np.sqrt(np.sum(df_err.loc[pos_mask , bin_cols].astype(float).values**2, axis=0))

        Err_mc = np.zeros(len(bin_cols))
        if np.nansum(df_err.loc[mc_mask, bin_cols]) > 0:    
            Err_mc   = df_err.loc[mc_mask, bin_cols].astype(float).values.flatten()

        # Positron Subtraction
        Yield_sub = Yield_elec - Yield_pos
        Err_sub   = np.sqrt(Err_elec**2 + Err_pos**2)

        # Ratio = (e− − e+) / MC
        ratio = np.divide(Yield_sub, Yield_mc, out=np.full_like(Yield_sub, np.nan), where=(Yield_mc > 0))
        ratio_err = np.full_like(ratio, np.nan)

        good = (Yield_sub > 0) & (Yield_mc > 0)
        ratio_err[good] = ratio[good] * np.sqrt((Err_sub[good] / Yield_sub[good])**2 + (Err_mc[good] / Yield_mc[good])**2)

        # Bin centers
        bin_left_edges = np.array([float(cols) for cols in bin_cols])
        bin_width = bin_left_edges[1] - bin_left_edges[0]
        bin_centers = bin_left_edges + 0.5 * bin_width

        # Saving binned yield ratios to output CSV
        output_dir = f"{selected_target_shortname.upper()}"
        output_filepath = f"{output_dir}/DATA_to_MC_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}.csv"
        os.makedirs(output_dir, exist_ok=True)

        df_output = pd.DataFrame({
            "bin_center": np.round(bin_centers, 4),
            "Yield_sub": Yield_sub,
            "Err_sub": Err_sub,
            "Yield_mc": Yield_mc,
            "Err_mc": Err_mc,
            "ratio": ratio,
            "ratio_err": ratio_err
        })

        df_output.to_csv(output_filepath, index=False)
        print(f"Saved binned yield CSV → {output_filepath}")

        # -------------------------
        # Plotting
        # -------------------------
        valid = (np.isfinite(Yield_sub) & np.isfinite(Yield_mc) & np.isfinite(ratio))
        
        fig, (ax_top, ax_bot) = plt.subplots(2,1,figsize=(8.5,11/2), gridspec_kw={'height_ratios':[3,1]}, sharex=True)

        # Always plot main data and MC scatter points
        ax_top.errorbar(bin_centers, Yield_sub, yerr=Err_sub, fmt='o', markersize=3, color='navy', label="Data")
        ax_top.errorbar(bin_centers, Yield_mc,  yerr=Err_mc, fmt='o', markersize=3, color='red',  label="MC")

        # Enhanced info for H_gtr_dp only
        if var == "H_gtr_dp":
            Yield_elec_tot = np.nansum(Yield_elec)
            Yield_pos_tot  = np.nansum(Yield_pos)
            Yield_sub_tot  = np.nansum(Yield_sub)
            pos_frac   = 100.0 * Yield_pos_tot / Yield_elec_tot if Yield_elec_tot > 0 else np.nan
            
            ax_top.errorbar(bin_centers[valid], Yield_elec[valid], yerr=Err_elec[valid], fmt='o', markersize=3, color='darkorange', label='e⁻ yield')
            ax_top.errorbar(bin_centers[valid], Yield_pos[valid],  yerr=Err_pos[valid],  fmt='o', markersize=3, color='darkgreen',  label='e⁺ yield')
            
            textstr = (f"Total e⁻ yield: {Yield_elec_tot:8.2f}\n"
                       f"Total e⁺ yield: {Yield_pos_tot:8.2f}\n"
                       f"Data (e⁻ - e⁺) : {Yield_sub_tot:8.2f}\n"
                       f"Total MC: {np.nansum(Yield_mc):8.2f}\n"
                       f"e⁺ fraction: {pos_frac:6.2f} %\n"
                       f"Data/MC ratio: {100.0 * np.nansum(Yield_sub)/np.nansum(Yield_mc):6.2f} %")

            ax_top.text(0.02, 0.75, textstr, transform=ax_top.transAxes, fontsize=10, verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", alpha=0.9))
            
            plt.tight_layout()

        # Bottom panel: ratio
        ax_bot.errorbar(bin_centers[valid], ratio[valid], yerr=ratio_err[valid], fmt='o', markersize=3, color='darkmagenta')
        ax_bot.axhline(1, color='gray', linestyle='--')
        ax_bot.set_ylabel("Data/MC"); ax_bot.set_ylim(0.75,1.25)
        ax_bot.set_xlabel(var); ax_bot.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        ax_bot.grid()

        ax_top.set_ylabel("Charge-Normalized Weighted Counts")
        ax_top.set_title(f"{selected_run_type.upper()} {selected_beam_pass}Pass {selected_target_titlename} {var}: Data (e⁻ - e⁺) to MC Yields")
        ax_top.grid(); ax_top.legend(loc='upper right')

        fig.tight_layout(); fig.subplots_adjust(top=0.93, bottom=0.12, hspace=0.05)
        pp.savefig(fig); plt.close(fig)

    pp.close()
    print(f"\nSaved PDF → {pdf_output}")
    exit(1)    
    
# -----------------------------------------------------
# TESTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Loop over variables and produce figures for LIQUID targets
# -----------------------------------------------------
if selected_target_shortname in {"ld2", "lh2"}:
    dummy_dir = f"../MAKE_csvs/DUMMY"
    R_dummy = 1 / 3.7825 if selected_target_shortname == 'ld2' else 1 / 3.5274
    ordered_variables = ["H_gtr_dp"] + [v for v in variables if v != "H_gtr_dp"]
    
    for var in ordered_variables:
        histo_filepath = f"{base_dir}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_histo.csv"
        err_filepath   = f"{base_dir}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_err.csv"
        dummy_histo_filepath = f"{dummy_dir}/{selected_run_type}_{selected_beam_pass}pass_dummy_{var}_histo.csv"
        dummy_err_filepath   = f"{dummy_dir}/{selected_run_type}_{selected_beam_pass}pass_dummy_{var}_err.csv"

        df_yield = pd.read_csv(histo_filepath)
        df_err = pd.read_csv(err_filepath)
        df_yield_dummy = pd.read_csv(dummy_histo_filepath)
        df_err_dummy = pd.read_csv(dummy_err_filepath)
        
        runinfo_cols = ['runnum','charge','polarity']
        bin_cols = [c for c in df_yield.columns if c not in runinfo_cols]
        bin_cols_dummy = [c for c in df_yield_dummy.columns if c not in runinfo_cols]

        # Masks
        elec_mask = df_yield['polarity'] == '-'
        pos_mask = df_yield['polarity'] == '+'
        mc_mask = df_yield['polarity'] == 'mc'
        elec_mask_dummy = df_yield_dummy['polarity'] == '-'
        pos_mask_dummy = df_yield_dummy['polarity'] == '+'

        
        # Charge Normalization
        charge_norm_elec = 0.0
        if np.nansum(df_yield.loc[elec_mask, 'charge']) > 0:
            charge_norm_elec = 1000 / df_yield.loc[elec_mask, 'charge'].sum()
            df_yield.loc[elec_mask, bin_cols] *= charge_norm_elec
            df_err.loc[elec_mask, bin_cols] *= charge_norm_elec
        
        charge_norm_pos = 0.0
        if np.nansum(df_yield.loc[pos_mask, 'charge']) > 0:
            charge_norm_pos  = 1000 / df_yield.loc[pos_mask , 'charge'].sum()
            df_yield.loc[pos_mask , bin_cols] *= charge_norm_pos
            df_err.loc[pos_mask , bin_cols] *= charge_norm_pos

        charge_norm_elec_dummy = 0.0
        if np.nansum(df_yield_dummy.loc[elec_mask_dummy, 'charge']) > 0:
            charge_norm_elec_dummy = 1000 / df_yield_dummy.loc[elec_mask_dummy, 'charge'].sum()
            df_yield_dummy.loc[elec_mask_dummy, bin_cols_dummy] *= charge_norm_elec_dummy
            df_err_dummy.loc[elec_mask_dummy, bin_cols_dummy] *= charge_norm_elec_dummy

        charge_norm_pos_dummy = 0.0
        if np.nansum(df_yield_dummy.loc[pos_mask_dummy, 'charge']) > 0:
            charge_norm_pos_dummy = 1000 / df_yield_dummy.loc[pos_mask_dummy, 'charge'].sum()
            df_yield_dummy.loc[pos_mask_dummy, bin_cols_dummy] *= charge_norm_pos_dummy
            df_err_dummy.loc[pos_mask_dummy, bin_cols_dummy] *= charge_norm_pos_dummy
        
        # Sum Normalized Yields
        Yield_elec = np.zeros(len(bin_cols))
        if np.nansum(df_yield.loc[elec_mask, bin_cols]) > 0:
            Yield_elec = df_yield.loc[elec_mask, bin_cols].astype(float).sum(axis=0).values
            
        Yield_pos = np.zeros(len(bin_cols))
        if np.nansum(df_yield.loc[pos_mask, bin_cols]) > 0:
            Yield_pos = df_yield.loc[pos_mask , bin_cols].astype(float).sum(axis=0).values

        Yield_mc = np.zeros(len(bin_cols))
        if np.nansum(df_yield.loc[mc_mask, bin_cols]) > 0:
            Yield_mc = df_yield.loc[mc_mask ,  bin_cols].astype(float).values.flatten()

        Yield_elec_dummy = np.zeros(len(bin_cols_dummy))
        if np.nansum(df_yield_dummy.loc[elec_mask_dummy, bin_cols_dummy]) > 0:
            Yield_elec_dummy = df_yield_dummy.loc[elec_mask_dummy, bin_cols_dummy].astype(float).sum(axis=0).values
            
        Yield_pos_dummy = np.zeros(len(bin_cols_dummy))
        if np.nansum(df_yield_dummy.loc[pos_mask_dummy, bin_cols_dummy]) > 0:
            Yield_pos_dummy = df_yield_dummy.loc[pos_mask_dummy, bin_cols_dummy].astype(float).sum(axis=0).values

        # Error Propagation
        Err_elec = np.zeros(len(bin_cols))
        if np.nansum(df_err.loc[elec_mask, bin_cols]) > 0:
            Err_elec = np.sqrt(np.sum(df_err.loc[elec_mask, bin_cols].astype(float).values**2, axis=0))
            
        Err_pos = np.zeros(len(bin_cols))
        if np.nansum(df_err.loc[pos_mask, bin_cols]) > 0:    
            Err_pos  = np.sqrt(np.sum(df_err.loc[pos_mask , bin_cols].astype(float).values**2, axis=0))

        Err_mc = np.zeros(len(bin_cols))
        if np.nansum(df_err.loc[mc_mask, bin_cols]) > 0:    
            Err_mc   = df_err.loc[mc_mask, bin_cols].astype(float).values.flatten()

        Err_elec_dummy = np.zeros(len(bin_cols_dummy))
        if np.nansum(df_err_dummy.loc[elec_mask_dummy, bin_cols_dummy]) > 0:
            Err_elec_dummy = np.sqrt(np.sum(df_err_dummy.loc[elec_mask_dummy, bin_cols_dummy].astype(float).values**2, axis=0))
            
        Err_pos_dummy = np.zeros(len(bin_cols_dummy))
        if np.nansum(df_err_dummy.loc[pos_mask_dummy, bin_cols_dummy]) > 0:    
            Err_pos_dummy  = np.sqrt(np.sum(df_err_dummy.loc[pos_mask_dummy , bin_cols_dummy].astype(float).values**2, axis=0))                      
            
        # Positron and Dummy Subtraction
        Yield_sub = Yield_elec - Yield_pos - R_dummy * (Yield_elec_dummy - Yield_pos_dummy)
        Err_sub = np.sqrt(Err_elec**2 + Err_pos**2 + (R_dummy * Err_elec_dummy)**2 + (R_dummy * Err_pos_dummy)**2)
        
        # Ratio = (e− − e+) / MC
        ratio = np.divide(Yield_sub, Yield_mc, out=np.full_like(Yield_sub, np.nan), where=(Yield_mc > 0))
        ratio_err = np.full_like(ratio, np.nan)

        good = (Yield_sub > 0) & (Yield_mc > 0)
        ratio_err[good] = ratio[good] * np.sqrt((Err_sub[good] / Yield_sub[good])**2 + (Err_mc[good] / Yield_mc[good])**2)

        # Bin centers
        bin_left_edges = np.array([float(cols) for cols in bin_cols])
        bin_width = bin_left_edges[1] - bin_left_edges[0]
        bin_centers = bin_left_edges + 0.5 * bin_width

        # Saving binned yield ratios to output CSV
        output_dir = f"{selected_target_shortname.upper()}"
        output_filepath = f"{output_dir}/DATA_to_MC_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}.csv"
        os.makedirs(output_dir, exist_ok=True)

        df_output = pd.DataFrame({
            "bin_center": np.round(bin_centers, 4),
            "Yield_sub": Yield_sub,
            "Err_sub": Err_sub,
            "Yield_mc": Yield_mc,
            "Err_mc": Err_mc,
            "ratio": ratio,
            "ratio_err": ratio_err
        })

        df_output.to_csv(output_filepath, index=False)
        print(f"Saved binned yield CSV → {output_filepath}")

        # -------------------------
        # Plotting
        # -------------------------
        valid = (np.isfinite(Yield_sub) & np.isfinite(Yield_mc) & np.isfinite(ratio))
        
        fig, (ax_top, ax_bot) = plt.subplots(2,1,figsize=(8.5,11/2), gridspec_kw={'height_ratios':[3,1]}, sharex=True)

        # Always plot main data and MC scatter points
        ax_top.errorbar(bin_centers, Yield_sub, yerr=Err_sub, fmt='o', markersize=3, color='navy', label="Data")
        ax_top.errorbar(bin_centers, Yield_mc,  yerr=Err_mc, fmt='o', markersize=3, color='red',  label="MC")

        # Enhanced info for H_gtr_dp only
        if var == "H_gtr_dp":
            Yield_elec_tot = np.nansum(Yield_elec)
            Yield_pos_tot  = np.nansum(Yield_pos)
            Yield_sub_tot  = np.nansum(Yield_sub)
            pos_frac   = 100.0 * Yield_pos_tot / Yield_elec_tot if Yield_elec_tot > 0 else np.nan
            
            ax_top.errorbar(bin_centers[valid], Yield_elec[valid], yerr=Err_elec[valid], fmt='o', markersize=3, color='darkorange', label='e⁻ yield')
            ax_top.errorbar(bin_centers[valid], Yield_pos[valid],  yerr=Err_pos[valid],  fmt='o', markersize=3, color='darkgreen',  label='e⁺ yield')
            
            textstr = (f"Total e⁻ yield: {Yield_elec_tot:8.2f}\n"
                       f"Total e⁺ yield: {Yield_pos_tot:8.2f}\n"
                       f"Data (e⁻ - e⁺) :     {Yield_sub_tot:8.2f}\n"
                       f"Total MC:       {np.nansum(Yield_mc):8.2f}\n"
                       f"e⁺ fraction:    {pos_frac:6.2f} %\n"
                       f"Data/MC ratio:  {100.0 * np.nansum(Yield_sub)/np.nansum(Yield_mc):6.2f} %")

            ax_top.text(0.02, 0.75, textstr, transform=ax_top.transAxes, fontsize=10,
                        verticalalignment='top', bbox=dict(boxstyle="round", facecolor="white", alpha=0.9))
            
            plt.tight_layout()

        # Bottom panel: ratio
        ax_bot.errorbar(bin_centers[valid], ratio[valid], yerr=ratio_err[valid], fmt='o', markersize=3, color='darkmagenta')
        ax_bot.axhline(1, color='gray', linestyle='--')
        ax_bot.set_ylabel("Data/MC"); ax_bot.set_ylim(0.75,1.25)
        ax_bot.set_xlabel(var); ax_bot.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        ax_bot.grid()

        ax_top.set_ylabel("Charge-Normalized Weighted Counts")
        ax_top.set_title(f"{selected_run_type.upper()} {selected_beam_pass}Pass {selected_target_titlename} {var}: Data (e⁻ - e⁺) to MC Yields")
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
