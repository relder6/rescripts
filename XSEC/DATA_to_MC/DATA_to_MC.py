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
selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, selected_target_titlename, selected_target_A, selected_target_Z = get_common_run_inputs()

# -----------------------------------------------------
# Files
# -----------------------------------------------------
base_dir = f"../MAKE_csvs/{selected_target_shortname.upper()}"
    
variables = ["H_gtr_dp", "H_gtr_ph", "H_gtr_th", "H_kin_Q2", "H_kin_x_bj", "H_kin_W", "H_gtr_p", "H_gtr_y", "H_gtr_th", "H_gtr_ph", "H_dc_x_fp", "H_dc_xp_fp", "H_dc_y_fp", "H_dc_yp_fp"]

pdf_output = f"PDFs/DATA_to_MC_{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}.pdf"
pp = PdfPages(pdf_output)
    
# -----------------------------------------------------
# Loop over variables and produce figures for ALL targets
# -----------------------------------------------------
if selected_target_shortname in {"c", "cu", "al", "ld2", "lh2"}:
    ordered_variables = ["H_gtr_dp"] + [v for v in variables if v != "H_gtr_dp"]
    if selected_target_shortname in {"ld2", "lh2"}:
        dummy_dir = f"../MAKE_csvs/DUMMY"
    
    for var in ordered_variables:
        histo_filepath = f"{base_dir}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_histo.csv"
        err_filepath   = f"{base_dir}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_err.csv"
        if selected_target_shortname in {"ld2", "lh2"}:
            dummy_histo_filepath = f"{dummy_dir}/{selected_run_type}_{selected_beam_pass}pass_dummy_{var}_{selected_target_shortname}_histo.csv"
            dummy_err_filepath   = f"{dummy_dir}/{selected_run_type}_{selected_beam_pass}pass_dummy_{var}_{selected_target_shortname}_err.csv"

        df_yield = pd.read_csv(histo_filepath)
        df_err = pd.read_csv(err_filepath)
        if selected_target_shortname in {"ld2", "lh2"}:
            df_yield_dummy = pd.read_csv(dummy_histo_filepath)
            df_err_dummy = pd.read_csv(dummy_err_filepath)
        
        runinfo_cols = ['runnum','charge','polarity']
        bin_cols = [c for c in df_yield.columns if c not in runinfo_cols]
        if selected_target_shortname in {"ld2", "lh2"}:
            bin_cols_dummy = [c for c in df_yield_dummy.columns if c not in runinfo_cols]

        # Masks
        elec_mask = df_yield['polarity'] == '-'
        pos_mask = df_yield['polarity'] == '+'
        mc_mask = df_yield['polarity'] == 'mc'
        if selected_target_shortname in {"ld2", "lh2"}:
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
            
        if selected_target_shortname in {"ld2", "lh2"}:
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
            
        if selected_target_shortname in {"ld2", "lh2"}:
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

        if selected_target_shortname in {"ld2", "lh2"}:
            Err_elec_dummy = np.zeros(len(bin_cols_dummy))
            if np.nansum(df_err_dummy.loc[elec_mask_dummy, bin_cols_dummy]) > 0:
                Err_elec_dummy = np.sqrt(np.sum(df_err_dummy.loc[elec_mask_dummy, bin_cols_dummy].astype(float).values**2, axis=0))
            
            Err_pos_dummy = np.zeros(len(bin_cols_dummy))
            if np.nansum(df_err_dummy.loc[pos_mask_dummy, bin_cols_dummy]) > 0:    
                Err_pos_dummy  = np.sqrt(np.sum(df_err_dummy.loc[pos_mask_dummy , bin_cols_dummy].astype(float).values**2, axis=0))                      
            
        # Positron and Dummy Subtraction
        if selected_target_shortname in {"c", "cu", "al"}:
            Yield_sub = Yield_elec - Yield_pos
            Err_sub   = np.sqrt(Err_elec**2 + Err_pos**2)
        if selected_target_shortname in {"ld2", "lh2"}:
            Yield_sub = Yield_elec - Yield_pos - Yield_elec_dummy + Yield_pos_dummy
            Err_sub = np.sqrt(Err_elec**2 + Err_pos**2 + Err_elec_dummy**2 + Err_pos_dummy**2)
        
        # Ratio = (e− − e+) / MC
        ratio = np.divide(Yield_sub, Yield_mc, out=np.full_like(Yield_sub, np.nan), where=(Yield_mc > 0))
        ratio_err = np.full_like(ratio, np.nan)

        good = (Yield_sub > 0) & (Yield_mc > 0)
        ratio_err[good] = ratio[good] * np.sqrt((Err_sub[good] / Yield_sub[good])**2 + (Err_mc[good] / Yield_mc[good])**2)

        # Bin centers
        bin_centers = np.array([float(cols) for cols in bin_cols])

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

        # Enhanced info for H_gtr_dp only
        if var == "H_gtr_dp":
            Yield_tot = np.nansum(Yield_elec) + np.nansum(Yield_pos)
            Yield_elec_tot = np.nansum(Yield_elec)
            Yield_elec_frac = Yield_elec_tot / Yield_tot
            Yield_pos_tot  = np.nansum(Yield_pos)
            Yield_pos_frac = Yield_pos_tot / Yield_tot
            Yield_sub_tot  = np.nansum(Yield_sub)
            Yield_elec_plot = np.where(valid, Yield_elec, np.nan)
            Err_elec_plot   = np.where(valid, Err_elec, np.nan)
            Yield_pos_plot  = np.where(valid, Yield_pos, np.nan)
            Err_pos_plot    = np.where(valid, Err_pos, np.nan)

            ax_top.errorbar(bin_centers, Yield_sub, yerr=Err_sub, fmt='o', markersize=3, color='navy', label=f"Data (e⁻ - e⁺): {Yield_sub_tot:8.2f}")
            ax_top.errorbar(bin_centers, Yield_mc,  yerr=Err_mc, fmt='o', markersize=3, color='red',  label=f"MC: {np.nansum(Yield_mc):8.2f}")
            ax_top.errorbar(bin_centers, Yield_elec_plot, yerr=Err_elec_plot, fmt='o', markersize=3, color='darkorange', label=f"e⁻ yield: {Yield_elec_tot:8.2f}")
            ax_top.errorbar(bin_centers, Yield_pos_plot,  yerr=Err_pos_plot,  fmt='o', markersize=3, color='darkgreen',  label=f"e⁺ yield: {Yield_pos_tot:8.2f}")

            # Add a "dummy plot" that will only appear in the legend as a horizontal line
            ax_top.plot([0], [0], linestyle="", color="none", label=f"Total (e⁻ + e⁺): {np.nansum(Yield_elec)+np.nansum(Yield_pos):.2f}")
            ax_top.plot([0], [0], linestyle="", color="none", label=f"e⁺ / (e⁻ + e⁺): {100*Yield_pos_frac:.2f}%")
            ax_top.plot([0], [0], linestyle="", color="none", label=f"Data/MC ratio: {100*np.nansum(Yield_sub)/np.nansum(Yield_mc):.2f}%")

            handles, labels = ax_top.get_legend_handles_labels()
            order = [3, 4, 5, 6, 0, 1, 2]
            ax_top.legend([handles[i] for i in order], [labels[i] for i in order], loc='best', fontsize=10)

            # ax_bot.set_ylim(0.90,1.10)

        else:
            ax_top.errorbar(bin_centers, Yield_sub, yerr=Err_sub, fmt='o', markersize=3, color='navy', label="Data")
            ax_top.errorbar(bin_centers, Yield_mc,  yerr=Err_mc, fmt='o', markersize=3, color='red',  label="MC")
            ax_top.legend(loc='best', fontsize=10)

            ax_bot.set_ylim(0.75,1.25)

        # Bottom panel: ratio
        ax_bot.errorbar(bin_centers[valid], ratio[valid], yerr=ratio_err[valid], fmt='o', markersize=3, color='darkmagenta', label = 'Data/MC Ratio')
        ax_bot.axhline(1, color='gray', linestyle='--')
        ax_bot.set_ylabel("Data/MC")
        ax_bot.set_xlabel(var); ax_bot.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        ax_bot.grid()

        ax_top.set_ylabel("Charge-Normalized Weighted Counts")
        ax_top.set_title(f"{selected_run_type.upper()} {selected_beam_pass}Pass {selected_target_titlename} {var}: Data (e⁻ - e⁺) to MC")
        ax_top.grid()

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
