#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import os, sys, re
import mplhep

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import parse_run_type, parse_beam_pass, parse_target, parse_phase

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None
arg3 = sys.argv[3] if len(sys.argv) > 3 else None
arg4 = sys.argv[4] if len(sys.argv) > 4 else None

selected_run_type = parse_run_type(arg1)
selected_beam_pass, _ = parse_beam_pass(arg2)
target_abbrev, target_longname, target_shortname, target_A, target_Z = parse_target(arg3)
phase = parse_phase(arg4)

# -----------------------------------------------------
# Files
# -----------------------------------------------------
base_dir = f"../MAKE_csvs/{target_shortname.upper()}"
    
variables = ["H_gtr_dp", "H_gtr_ph", "H_gtr_th", "H_kin_Q2", "H_kin_x_bj", "H_kin_W", "H_gtr_p", "H_gtr_y", "H_dc_x_fp", "H_dc_xp_fp", "H_dc_y_fp", "H_dc_yp_fp", "H_kin_W2"]

pdf_dir = f"PDFs"

pdf_output = f"{pdf_dir}/DATA_to_MC_{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}.pdf"

ordered_variables = ["H_gtr_dp"] + [v for v in variables if v != "H_gtr_dp"]

input_filepath = (f"{base_dir}/{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}.csv")

if not os.path.isfile(input_filepath):
    print("****************************************************")
    print(f"No input CSV file found, {selected_run_type} {selected_beam_pass}pass {target_shortname} phase{phase}; skipping PDF creation.\n")
    sys.exit(0)

if target_abbrev in {"ld2", "lh2"}:
    dummy_dir = "../MAKE_csvs/DUMMY"
    dummy_filepath = f"{dummy_dir}/{selected_run_type}_{selected_beam_pass}pass_phase{phase}_dummy.csv"

    if not os.path.isfile(dummy_filepath):
        print("****************************************************")
        print(f"No dummy CSV file found, {selected_run_type} {selected_beam_pass}pass dummy phase{phase}; skipping PDF creation.\n")
        sys.exit(0)
        
pp = PdfPages(pdf_output)
    
# -----------------------------------------------------
# Loop over variables and produce figures for ALL targets
# -----------------------------------------------------
if target_abbrev in {"al", "c", "cu", "ld2", "lh2", "dummy_up", "dummy_down"}:

    df_all = pd.read_csv(input_filepath)

    for var in ordered_variables:
        df_var = df_all[df_all["variable"] == var].copy()

        if df_var.empty:
            continue

        nbins = int(df_var.iloc[0]["nbins"])
        bin_cols = [f"bin{i}" for i in range(nbins)]
        
        if not bin_cols:
            continue
        
        
        if target_abbrev in {"ld2", "lh2"}:
            dummy_filepath = f"{dummy_dir}/{selected_run_type}_{selected_beam_pass}pass_phase{phase}_dummy.csv"
            df_dummy = pd.read_csv(dummy_filepath)
            df_var_dummy = df_dummy[df_dummy["variable"] == var]

            if df_var_dummy.empty:
                continue
            
            nbins_dummy = int(df_var_dummy.iloc[0]["nbins"])

            bin_cols_dummy = [f"bin{i}" for i in range(nbins_dummy)]

            if not bin_cols_dummy:
                continue
        
        runinfo_cols = ["variable","type","runnum","charge","polarity","nbins","bin_min","bin_max"]

        # Masks
        elec_mask = (df_var["type"] == "data") & (df_var["polarity"] == "-")
        pos_mask  = (df_var["type"] == "data") & (df_var["polarity"] == "+")
        elec_mask_err = (df_var["type"] == "err") & (df_var["polarity"] == "-")
        pos_mask_err  = (df_var["type"] == "err") & (df_var["polarity"] == "+")
        
        mc_mask = (df_var["type"] == "mc") & (df_var["polarity"] == "-")
        mc_mask_err = (df_var["type"] == "mc_err") & (df_var["polarity"] == "-")
        
        if target_abbrev in {"ld2", "lh2"}:
            elec_mask_dummy = (df_var_dummy["type"] == "data") & (df_var_dummy["polarity"] == "-")
            pos_mask_dummy = (df_var_dummy["type"] == "data") & (df_var_dummy["polarity"] == "+")
            elec_mask_dummy_err = (df_var_dummy["type"] == "err") & (df_var_dummy["polarity"] == "-")
            pos_mask_dummy_err = (df_var_dummy["type"] == "err") & (df_var_dummy["polarity"] == "+")

        # Charge Normalization
        charge_norm_neg = 0.0
        if np.nansum(df_var.loc[elec_mask, "charge"]) > 0:
            charge_norm_neg = 1000 / df_var.loc[elec_mask, "charge"].sum()
            df_var.loc[elec_mask, bin_cols] *= charge_norm_neg
        
        charge_norm_pos = 0.0
        if np.nansum(df_var.loc[pos_mask, "charge"]) > 0:
            charge_norm_pos  = 1000 / df_var.loc[pos_mask , "charge"].sum()
            df_var.loc[pos_mask , bin_cols] *= charge_norm_pos
            
        if target_abbrev in {"ld2", "lh2"}:
            charge_norm_neg_dummy = 0.0
            if np.nansum(df_var_dummy.loc[elec_mask_dummy, "charge"]) > 0:
                charge_norm_neg_dummy = 1000 / df_var_dummy.loc[elec_mask_dummy, "charge"].sum()
                df_var_dummy.loc[elec_mask_dummy, bin_cols_dummy] *= charge_norm_neg_dummy

            charge_norm_pos_dummy = 0.0
            if np.nansum(df_var_dummy.loc[pos_mask_dummy, "charge"]) > 0:
                charge_norm_pos_dummy = 1000 / df_var_dummy.loc[pos_mask_dummy, "charge"].sum()
                df_var_dummy.loc[pos_mask_dummy, bin_cols_dummy] *= charge_norm_pos_dummy
        
        # Sum Normalized Yields
        yield_neg = np.zeros(len(bin_cols))
        if np.nansum(df_var.loc[elec_mask, bin_cols]) > 0:
            yield_neg = df_var.loc[elec_mask, bin_cols].astype(float).sum(axis=0).values
            
        yield_pos = np.zeros(len(bin_cols))
        if np.nansum(df_var.loc[pos_mask, bin_cols]) > 0:
            yield_pos = df_var.loc[pos_mask , bin_cols].astype(float).sum(axis=0).values

        yield_mc = np.zeros(len(bin_cols))
        if np.nansum(df_var.loc[mc_mask, bin_cols]) > 0:
            yield_mc = df_var.loc[mc_mask, bin_cols].astype(float).sum(axis=0).values
            
        if target_abbrev in {"ld2", "lh2"}:
            yield_neg_dummy = np.zeros(len(bin_cols_dummy))
            if np.nansum(df_var_dummy.loc[elec_mask_dummy, bin_cols_dummy]) > 0:
                yield_neg_dummy = df_var_dummy.loc[elec_mask_dummy, bin_cols_dummy].astype(float).sum(axis=0).values
            
            yield_pos_dummy = np.zeros(len(bin_cols_dummy))
            if np.nansum(df_var_dummy.loc[pos_mask_dummy, bin_cols_dummy]) > 0:
                yield_pos_dummy = df_var_dummy.loc[pos_mask_dummy, bin_cols_dummy].astype(float).sum(axis=0).values

        # Error Propagation
        err_neg = df_var.loc[elec_mask_err, bin_cols].astype(float).iloc[0].values
        err_pos = df_var.loc[pos_mask_err,  bin_cols].astype(float).iloc[0].values
        err_mc  = df_var.loc[mc_mask_err,   bin_cols].astype(float).iloc[0].values

        err_neg *= charge_norm_neg
        err_pos *= charge_norm_pos
        # err_neg = np.zeros(len(bin_cols))
        # if np.nansum(df_var.loc[elec_mask_err, bin_cols]) > 0:
        #     err_neg = np.sqrt(np.sum(df_var.loc[elec_mask_err, bin_cols].astype(float).values**2, axis=0))
            
        # err_pos = np.zeros(len(bin_cols))
        # if np.nansum(df_var.loc[pos_mask_err, bin_cols]) > 0:    
        #     err_pos  = np.sqrt(np.sum(df_var.loc[pos_mask_err, bin_cols].astype(float).values**2, axis=0))

        # err_mc = np.zeros(len(bin_cols))
        # if np.nansum(df_var.loc[mc_mask_err, bin_cols]) > 0:    
        #     err_mc = np.sqrt(np.sum(df_var.loc[mc_mask_err, bin_cols].astype(float).values**2, axis=0))

        if target_abbrev in {"ld2", "lh2"}:
            err_neg_dummy = df_var_dummy.loc[elec_mask_dummy_err, bin_cols_dummy].astype(float).iloc[0].values
            err_pos_dummy = df_var_dummy.loc[pos_mask_dummy_err, bin_cols_dummy].astype(float).iloc[0].values

            err_neg_dummy *= charge_norm_neg_dummy
            err_pos_dummy *= charge_norm_pos_dummy
            
            # err_neg_dummy = np.zeros(len(bin_cols_dummy))
            # if np.nansum(df_var_dummy.loc[elec_mask_dummy_err, bin_cols_dummy]) > 0:
            #     err_neg_dummy = np.sqrt(np.sum(df_var_dummy.loc[elec_mask_dummy_err, bin_cols_dummy].astype(float).values**2, axis=0))
            
            # err_pos_dummy = np.zeros(len(bin_cols_dummy))
            # if np.nansum(df_var_dummy.loc[pos_mask_dummy_err, bin_cols_dummy]) > 0:    
            #     err_pos_dummy  = np.sqrt(np.sum(df_var_dummy.loc[pos_mask_dummy_err, bin_cols_dummy].astype(float).values**2, axis=0))     
                
            
        # Positron and Dummy Subtraction
        if target_abbrev in {"al", "c", "cu", "dummy_up", "dummy_down"}:
            yield_sub = yield_neg - yield_pos
            err_sub   = np.sqrt(err_neg**2 + err_pos**2)
        if target_abbrev in {"ld2", "lh2"}:
            yield_sub = yield_neg - yield_pos - yield_neg_dummy + yield_pos_dummy
            err_sub = np.sqrt(err_neg**2 + err_pos**2 + err_neg_dummy**2 + err_pos_dummy**2)
        
        # Ratio = (e− − e+) / MC
        ratio = np.divide(yield_sub, yield_mc, out=np.full_like(yield_sub, np.nan), where=(yield_mc > 0))
        ratio_err = np.full_like(ratio, np.nan)

        good = (yield_sub > 0) & (yield_mc > 0)
        ratio_err[good] = ratio[good] * np.sqrt((err_sub[good] / yield_sub[good])**2 + (err_mc[good] / yield_mc[good])**2)

        # Bin centers
        # bin_centers = np.array([float(cols) for cols in bin_cols])
        nbins = int(df_var.iloc[0]["nbins"])
        bin_min = df_var.iloc[0]["bin_min"]
        bin_max = df_var.iloc[0]["bin_max"]

        bin_centers = np.linspace(bin_min,bin_max,nbins,endpoint=False) + (bin_max-bin_min)/(2*nbins)

        # Saving binned yield ratios to output CSV
        output_dir = f"{target_shortname.upper()}"
        output_filepath = f"{output_dir}/DATA_to_MC_{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}_{var}.csv"
        os.makedirs(output_dir, exist_ok=True)

        df_output = pd.DataFrame({"bin_center": np.round(bin_centers, 4),
                                  "yield_sub": yield_sub,
                                  "err_sub": err_sub,
                                  "yield_mc": yield_mc,
                                  "err_mc": err_mc,
                                  "ratio": ratio,
                                  "ratio_err": ratio_err})

        df_output.to_csv(output_filepath, index=False)

        # -------------------------
        # Plotting
        # -------------------------
        valid = (np.isfinite(yield_sub) & np.isfinite(yield_mc) & np.isfinite(ratio))

        plt.style.use(mplhep.style.ROOT)
        plt.rcParams.update({"figure.titlesize": 14,
                             "axes.titlesize": 12,
                             "axes.labelsize": 10,
                             "legend.fontsize": 8,
                             "xtick.labelsize": 8,
                             "ytick.labelsize": 8})
        
        # fig, (ax_top, ax_bot) = plt.subplots(2,1,figsize=(8.5, 6), gridspec_kw={"height_ratios":[3,1]}, sharex=True)
        fig, (ax_top, ax_bot) = plt.subplots(2, 1,figsize=(6, 5),gridspec_kw={"height_ratios":[3,1]},sharex=True,constrained_layout=True)
        # Enhanced info for H_gtr_dp only
        if var == "H_gtr_dp" or (var == "H_gtr_y" and target_abbrev in {"ld2", "lh2"}):
            yield_neg_tot = np.nansum(yield_neg) # Negative polarity run, this is raw counts including background
            yield_sub_tot  = np.nansum(yield_sub)
            yield_pos_tot  = np.nansum(yield_pos) # Positive polarity run, this is the background that should be subtracted
            yield_mc_tot = np.nansum(yield_mc)
            ratio_pos_to_neg = yield_pos_tot / yield_neg_tot
            ratio_data_to_mc = yield_sub_tot / yield_mc_tot

            err_neg_tot = np.sqrt(np.nansum(err_neg**2))
            err_pos_tot = np.sqrt(np.nansum(err_pos**2))
            err_sub_tot = np.sqrt(np.nansum(err_sub**2))
            err_mc_tot = np.sqrt(np.nansum(err_mc**2))
            err_ratio_pos_to_neg = ratio_pos_to_neg * np.sqrt((err_pos_tot / yield_pos_tot)**2 +(err_neg_tot / yield_neg_tot)**2)
            err_ratio_data_to_mc = ratio_data_to_mc * np.sqrt((err_sub_tot / yield_sub_tot)**2 + (err_mc_tot / yield_mc_tot)**2)

            if target_abbrev in {"ld2", "lh2"}:
                yield_dummy_tot = np.nansum(yield_neg_dummy - yield_pos_dummy)
                err_dummy_tot = np.sqrt(np.nansum(err_neg_dummy**2 + err_pos_dummy**2))
                dummy_plot = np.where(valid, yield_neg_dummy - yield_pos_dummy, np.nan)
                err_dummy_plot = np.where(valid, np.sqrt(err_neg_dummy**2 + err_pos_dummy**2), np.nan)
                data_label = f"$Yield_{{Sub}}$ ($e^- - e^-_{{cs}}$ - Dummy): {yield_sub_tot:8.2f} $\pm$ {err_sub_tot:.2f}"
            else:
                data_label =f"$Yield_{{Sub}}$ ($e^- - e^-_{{cs}})$: {yield_sub_tot:8.2f} $\pm$ {err_sub_tot:.2f}"
            
            yield_neg_plot = np.where(valid, yield_neg, np.nan)
            err_neg_plot   = np.where(valid, err_neg, np.nan)
            yield_pos_plot  = np.where(valid, yield_pos, np.nan)
            err_pos_plot    = np.where(valid, err_pos, np.nan)
            
            ax_top.errorbar(bin_centers, yield_neg_plot, yerr=err_neg_plot, fmt="o", markersize=3, color="darkorange",
                            label=f"$Yield_{{Raw}}$ ($e^- + e^-_{{cs}}$): {yield_neg_tot:8.2f} $\pm$ {err_neg_tot:.2f}")
            ax_top.errorbar(bin_centers, yield_sub, yerr=err_sub, fmt="o", markersize=3, color="navy", label=data_label)
            ax_top.errorbar(bin_centers, yield_mc,  yerr=err_mc, fmt="o", markersize=3, color="red",
                            label=f"$Yield_{{MC}}$: {yield_mc_tot:8.2f} $\pm$ {err_mc_tot:.2f}")
            ax_top.errorbar(bin_centers, yield_pos_plot,  yerr=err_pos_plot,  fmt="o", markersize=3, color="darkgreen",
                            label=f"$Yield_{{Pos.}} (e^+_{{cs}})$: {yield_pos_tot:8.2f} $\pm$ {err_pos_tot:.2f}")
            if target_abbrev in {"ld2", "lh2"}:
                ax_top.errorbar(bin_centers, dummy_plot, yerr=err_dummy_plot, fmt="o", markersize = 3, color = "lightseagreen",
                                label = f"$Yield_{{Dummy}}$: {yield_dummy_tot:8.2f} $\pm$ {err_dummy_tot:8.2f}")
            ax_top.plot([0], [0], linestyle="", color="none",
                        label=f"$Y_{{Sub}}/Y_{{MC}}$ ratio: {100*ratio_data_to_mc:8.2f}% $\pm$ {100*err_ratio_data_to_mc:.2f}%")
            ax_top.plot([0], [0], linestyle="", color="none",
                        label=f"$Y_{{Pos.}}/Y_{{Raw}}$: {100*ratio_pos_to_neg:.2f}% $\pm$ {100*err_ratio_pos_to_neg:.2f}%")

            handles_y, labels_y = ax_top.get_legend_handles_labels()

            ax_top.legend(handles=[*handles_y],loc="center left",frameon=True,facecolor="white",edgecolor="gray",framealpha=0.6,fontsize=9,handlelength=0)

            if var == "H_gtr_dp":
                print(f"\n****************************************************")
                print(f"{selected_run_type.upper()}\t{selected_beam_pass}Pass\t{target_longname}")
                print(f"\nYield_Raw:\t{yield_neg_tot:10.2f}\t±  {err_neg_tot:8.2f}")
                print(f"Yield_Sub:\t{yield_sub_tot:10.2f}\t±  {err_sub_tot:8.2f}")
                print(f"Yield_MC:\t{yield_mc_tot:10.2f}\t±  {err_mc_tot:8.2f}")
                print(f"Yield_Pos:\t{yield_pos_tot:10.2f}\t±  {err_pos_tot:8.2f}")
                print(f"\nY_Sub/Y_MC:\t{100*ratio_data_to_mc:8.2f}%\t±  {100*err_ratio_data_to_mc:8.2f}%")
                print(f"Y_Pos/Y_Raw:\t{100*ratio_pos_to_neg:8.2f}%\t±  {100*err_ratio_pos_to_neg:8.2f}%")
        else:
            ax_top.errorbar(bin_centers, yield_sub, yerr=err_sub, fmt="o", markersize=3, color="navy", label="$Yield_{{Sub}}$")
            ax_top.errorbar(bin_centers, yield_mc,  yerr=err_mc, fmt="o", markersize=3, color="red",  label="$Yield_{{MC}}$")
            ax_top.legend(loc="best", frameon = True, facecolor = "white", edgecolor = "gray", framealpha = 0.6, fontsize=9, handlelength = 0)

            ax_bot.set_ylim(0.75,1.25)
            
        ax_top.set_ylabel("Charge-Normalized Yields")
        ax_top.set_title(f"{selected_run_type.upper()} {selected_beam_pass}Pass {target_longname} {var}: $Yield_{{Sub}}$ to $Yield_{{MC}}$")
        ax_top.grid()
        
        ax_bot.errorbar(bin_centers[valid], ratio[valid], yerr=ratio_err[valid], fmt="o", markersize=3, color="darkmagenta", label = "Data/MC Ratio")
        # ax_bot.axhline(1, color="gray", linestyle="--")
        ax_bot.set_ylabel("$Y_{{Sub}}/Y_{{MC}}$")
        ax_bot.set_xlabel(var); ax_bot.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        ax_bot.grid()
        


        pp.savefig(fig)
        plt.close(fig)

    pp.close()
    print(f"\nSaved CSVs → {output_dir}/ and PDF → {pdf_dir}/")
    
    sys.exit(0)

# -----------------------------------------------------
# Exception for other types (HEE, HEEP, Optics, HOLE)
# -----------------------------------------------------
else:
    print(f"\nTarget {target_abbrev} not currently supported by this script.")
    exit(1)
