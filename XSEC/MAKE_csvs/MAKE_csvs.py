#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import boost_histogram as bh
import os, re, sys
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import get_common_run_inputs, get_data_cuts

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
USING_DELTA_CORR = True

selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, selected_target_titlename, selected_target_A, selected_target_Z = get_common_run_inputs()

# -----------------------------------------------------
# File paths
# -----------------------------------------------------
rootfile_dir = "/work/hallc/c-rsidis/skimfiles/pass0"
mc_dir = "/work/hallc/c-rsidis/relder/mc-single-arm"
input_settings_filepath = f"../../FILTER_type/{selected_target_shortname.upper()}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_runs.dat"

output_dir = f"{selected_target_shortname.upper()}"

# -----------------------------------------------------
# Reading in monte-carlo report to obtain normfac
# -----------------------------------------------------
if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
    mc_filepath = f"{mc_dir}/worksim/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}.root"
    if not os.path.exists(mc_filepath):
        print(f"WARNING:\tNo mc-single-arm generated root file found for {selected_run_type}_{selected_beam_pass}_{selected_target_shortname}. Exiting...")
        exit(1)
    mc_report_filepath = f"{mc_dir}/outfiles/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}.out"
    with open(mc_report_filepath, "r") as infile:
        normfac_line = [line for line in infile if "NORMFAC" in line.upper()]
    normfac = float(normfac_line[0].split(":")[1].split()[0]) if normfac_line else None

# -----------------------------------------------------
# Defining branches, using uproot to put them in data frames
# -----------------------------------------------------
runnums = []
weight = []
charge = []
polarity = []

with open(input_settings_filepath, "r") as infile:
    next(infile) # Skipping the header line here
    for line in infile:
        parts = line.strip().split("\t")
        runnums.append(parts[0])
        charge.append(parts[12])
        weight.append(parts[19])
        polarity.append(parts[6].strip()[0])

# -----------------------------------------------------
# Defining branches, using uproot to put them in data frames
# -----------------------------------------------------
branches = ["H_gtr_dp", "H_cal_etottracknorm", "H_gtr_ph",
            "H_gtr_th", "H_gtr_x", "H_gtr_y",
            "H_kin_Q2", "H_kin_x_bj", "H_kin_W",
            "H_cer_npeSum", "H_gtr_p",
            "H_gtr_y", "H_gtr_th", "H_gtr_ph",
            "H_dc_x_fp", "H_dc_xp_fp", "H_dc_y_fp", "H_dc_yp_fp"]

branches_mc = ["hsdelta", "hsyptar", "hsxptar", "q2", "xb", "w", "weight", "eprime", "hsytar", "hsxptar", "hsyptar", "hsxfp", "hsxpfp", "hsyfp", "hsypfp"]

variable_mc_map = {
    "H_gtr_dp": "hsdelta",
    "H_gtr_ph": "hsyptar",
    "H_gtr_th": "hsxptar",
    "H_kin_Q2": "q2",
    "H_kin_x_bj": "xb",
    "H_kin_W": "w",
    "H_gtr_p": "eprime",
    "H_gtr_y": "hsytar",
    "H_gtr_th": "hsxptar",
    "H_gtr_ph": "hsyptar",
    "H_dc_x_fp": "hsxfp",
    "H_dc_xp_fp": "hsxpfp",
    "H_dc_y_fp": "hsyfp",
    "H_dc_yp_fp": "hsypfp"
    }

# -----------------------------------------------------
# Binning
# -----------------------------------------------------
if selected_beam_pass == "4":
    custom_bins = {
        "H_gtr_dp": dict(binnum = 20, min = -10.000, max = 10.000),
        "H_gtr_ph": dict(binnum = 20, min = -0.050, max = 0.050),
        "H_gtr_th": dict(binnum = 20, min = -0.100, max = 0.100),
        "H_kin_Q2": dict(binnum = 20, min = 2.400, max = 4.200),
        "H_kin_x_bj": dict(binnum = 20, min = 0.175, max = 0.325),
        "H_kin_W": dict(binnum = 20, min = 3.100, max = 3.500),
        "H_gtr_p": dict(binnum = 100, min = 1.3, max = 1.8),
        "H_gtr_y": dict(binnum = 100, min = -4, max = 4),
        "H_gtr_th": dict(binnum = 100, min = -0.1, max = 0.1),
        "H_gtr_ph": dict(binnum = 100, min = -0.05, max = 0.05),
        "H_dc_x_fp": dict(binnum = 20, min = -50, max = 50),
        "H_dc_xp_fp": dict(binnum = 20, min = -0.08, max = 0.08),
        "H_dc_y_fp": dict(binnum = 20, min = -30, max = 30),
        "H_dc_yp_fp": dict(binnum = 20, min = -0.04, max = 0.04),
    }
    
if selected_beam_pass == "5":
    custom_bins = {
        "H_gtr_dp": dict(binnum = 20, min = -10.000, max = 10.000),
        "H_gtr_ph": dict(binnum = 20, min = -0.050, max = 0.050),
        "H_gtr_th": dict(binnum = 20, min = -0.100, max = 0.100),
        "H_kin_Q2": dict(binnum = 20, min = 2.400, max = 4.200),
        "H_kin_x_bj": dict(binnum = 20, min = 0.175, max = 0.325),
        "H_kin_W": dict(binnum = 20, min = 3.100, max = 3.500),
        "H_gtr_p": dict(binnum = 100, min = 3.2, max = 5.2),
        "H_gtr_y": dict(binnum = 100, min = -4, max = 4),
        "H_gtr_th": dict(binnum = 100, min = -0.1, max = 0.1),
        "H_gtr_ph": dict(binnum = 100, min = -0.05, max = 0.05),
        "H_dc_x_fp": dict(binnum = 20, min = -50, max = 50),
        "H_dc_xp_fp": dict(binnum = 20, min = -0.08, max = 0.08),
        "H_dc_y_fp": dict(binnum = 20, min = -30, max = 30),
        "H_dc_yp_fp": dict(binnum = 20, min = -0.04, max = 0.04),
    }

# -----------------------------------------------------
# Data histogram and csv creation
# -----------------------------------------------------
hist_data = {}
hist_err_data = {}
bin_edges_dict = {}

cuts = get_data_cuts()

for var, bins in custom_bins.items():
    axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"], underflow=True, overflow=True)
    bin_edges_dict[var] = axis.edges
    hist_data[var] = []
    hist_err_data[var] = []

for i, runnum in enumerate(runnums):
     rootfile_path = f"{rootfile_dir}/skimmed_hms_coin_replay_production_{runnum}_-1.root"
     if not os.path.exists(rootfile_path):
         print(f"WARNING: Missing {rootfile_path}, skipping...")
         continue
     df = pd.DataFrame(uproot.open(rootfile_path)["T"].arrays(branches, library = "np"))
     data_cut = (df["H_gtr_dp"].between(cuts["H_gtr_dp_min_cut"], cuts["H_gtr_dp_max_cut"]) &
                 (df["H_cer_npeSum"] > cuts["H_cer_npeSum_cut"]) &
                 (df["H_cal_etottracknorm"] > cuts["H_cal_etottracknorm_cut"]))
     df_cut = df[data_cut].copy()

     run_weight = float(weight[i])

     for var, bins in custom_bins.items():
         axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"], underflow=True, overflow=True)

         if selected_target_shortname == "dummy":

             hist_lh2 = bh.Histogram(axis, storage=bh.storage.Weight())
             hist_ld2 = bh.Histogram(axis, storage=bh.storage.Weight())

             yvals = df_cut["H_gtr_y"].values
             varvals = df_cut[var].values

             y_mid = 0.0
             upstream_mask = yvals >= y_mid
             downstream_mask = yvals < y_mid

             # R_dummy_lh2_up, R_dummy_lh2_down = 1 / 3.5274, 1 / 3.5274 # Averaging two walls
             # R_dummy_ld2_up, R_dummy_ld2_down = 1 / 3.7825, 1 / 3.7825 # Averaging two walls
             
             R_dummy_lh2_up, R_dummy_lh2_down = 1 / 4.04033, 1 / 3.12459
             R_dummy_ld2_up, R_dummy_ld2_down = 1 / 4.66192, 1 / 3.17445

             weights_lh2 = np.zeros(len(yvals))
             weights_lh2[upstream_mask] = R_dummy_lh2_up * run_weight
             weights_lh2[downstream_mask] = R_dummy_lh2_down * run_weight

             weights_ld2 = np.zeros(len(yvals))
             weights_ld2[upstream_mask] = R_dummy_ld2_up * run_weight
             weights_ld2[downstream_mask] = R_dummy_ld2_down * run_weight

             hist_lh2.fill(df_cut[var], weight=weights_lh2)
             hist_ld2.fill(df_cut[var], weight=weights_ld2)

             counts_lh2 = hist_lh2.view().value
             counts_ld2 = hist_ld2.view().value

             errors_lh2 = np.sqrt(hist_lh2.view().variance)
             errors_ld2 = np.sqrt(hist_ld2.view().variance)

             hist_data.setdefault(var + "_lh2", []).append([runnum, charge[i], polarity[i]] + counts_lh2.tolist())
             hist_err_data.setdefault(var + "_lh2", []).append([runnum, charge[i], polarity[i]] + errors_lh2.tolist())

             hist_data.setdefault(var + "_ld2", []).append([runnum, charge[i], polarity[i]] + counts_ld2.tolist())
             hist_err_data.setdefault(var + "_ld2", []).append([runnum, charge[i], polarity[i]] + errors_ld2.tolist())
             
         else:     
             hist = bh.Histogram(axis, storage=bh.storage.Weight())
             hist.fill(df_cut[var], weight=np.full(len(df_cut[var]), run_weight))

             counts = hist.view().value
             errors = np.sqrt(hist.view().variance)
         
             hist_data[var].append([runnum, charge[i], polarity[i]]+counts.tolist())
             hist_err_data[var].append([runnum, charge[i], polarity[i]]+errors.tolist())

# -----------------------------------------------------
# Monte carlo histogram and csv creation
# -----------------------------------------------------
if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
    mc_file = uproot.open(mc_filepath)
    mc_tree = mc_file["h10"]
    df_mc = pd.DataFrame(mc_tree.arrays(branches_mc, library="np"))
    mc_cut = (df_mc["hsdelta"].between(cuts["H_gtr_dp_min_cut"], cuts["H_gtr_dp_max_cut"]))
    df_mc_cut = df_mc[mc_cut].copy()
    # print("DEBUG: df_mc shape:", df_mc.shape)
    # print("DEBUG: df_mc columns:", df_mc.columns.tolist())
    # print("DEBUG: hsdelta stats:",
    #       "min=", np.min(df_mc["hsdelta"]),
    #       "max=", np.max(df_mc["hsdelta"]))

    # print("DEBUG: df_mc_cut shape:", df_mc_cut.shape)
    # if df_mc_cut.empty:
    #     print("ERROR: df_mc_cut is EMPTY — hsdelta cut removed everything.")
    mc_hist_data = {}
    mc_hist_err = {}

for var, bins in custom_bins.items():
    mc_var = None
    
    if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        mc_var = variable_mc_map.get(var)
        
    if mc_var is None or selected_target_shortname in {"dummy", "optics1", "optics2", "hole"}:
        continue
    
    axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"])
    if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        hist_mc = bh.Histogram(axis, storage=bh.storage.Weight())

        if "weight" in df_mc_cut.columns:
            deltatmp = df_mc_cut["hsdelta"].values

            h1 = 1.0069
            h2 = 0.34018e-02
            h3 = -0.71161e-03
            h4 = -0.12060e-04
            h5 = 0.11322e-04
            h6 = -0.78222e-06

            deltacorr = (h1 + h2 * deltatmp + h3 * deltatmp**2 + h4 * deltatmp**3 + h5 * deltatmp**4 + h6 * deltatmp**5)

            if USING_DELTA_CORR:
                event_weights = df_mc_cut["weight"].values * normfac * deltacorr
            else:
                event_weights = df_mc_cut["weight"].values * normfac
        else:
            print("Weights branch not found.  Exiting...")
            print(df_mc_cut.columns.tolist())
            exit(1)
        # print(f"\nDEBUG: var={var}, mc_var={mc_var}")
        # print("DEBUG: data min/max:", df_mc_cut[mc_var].min(), df_mc_cut[mc_var].max())
        # print("DEBUG: bins for this variable:", bins["min"], bins["max"])
        # print("DEBUG: weight min/max/sum:",event_weights.min(),event_weights.max(), event_weights.sum())
        hist_mc.fill(df_mc_cut[mc_var].values, weight=event_weights)
        # print("DEBUG: hist sum after fill:", hist_mc.sum())
        # print("DEBUG: first 10 bin contents:", hist_mc.view().value[:10])

        mc_hist_data[var] = hist_mc.view().value.tolist()
        mc_hist_err[var] = np.sqrt(hist_mc.view().variance).tolist()


for var, rows in hist_data.items():
    base_var = var.replace("_lh2", "").replace("_ld2", "")
    bin_edges = bin_edges_dict[base_var]
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_labels = [f"{i:.6f}" for i in bin_centers]
    columns = ["runnum", "charge", "polarity"] + bin_labels

    # Saving here for counts
    hist_df = pd.DataFrame(rows, columns=columns)
    if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        mc_row = ["MC", "0", "mc"] + mc_hist_data[var]
        hist_df.loc[0] = mc_row
        hist_df.index = range(len(hist_df))
    else:
        hist_df.index = range(len(hist_df))
    
    output_filename = f"{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_histo.csv"
    output_filepath = f"{output_dir}/{output_filename}"
    hist_df.to_csv(output_filepath, index=False)

    # Saving here for error bars
    hist_err_df = pd.DataFrame(hist_err_data[var], columns=columns)
    if selected_target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        mc_err_row = ["MC", "0", "0"] + mc_hist_err[var]
        hist_err_df.loc[0] = mc_err_row
   
    hist_err_df.index = range(len(hist_df))
    
    output_err_filename = f"{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_{var}_err.csv"
    output_err_filepath = f"{output_dir}/{output_err_filename}"
    hist_err_df.to_csv(output_err_filepath, index=False)
    
    print(f"Saved {output_filename} and {output_err_filename} to {output_filepath}.")
