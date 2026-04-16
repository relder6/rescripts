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
import csv
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import parse_run_type, parse_beam_pass, parse_target, get_data_cuts

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
USING_DELTA_CORR = True

arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None
arg3 = sys.argv[3] if len(sys.argv) > 3 else None

selected_run_type = parse_run_type(arg1)
selected_beam_pass, beam_prefix = parse_beam_pass(arg2)
target_abbrev, target_longname, target_shortname, target_A, target_Z = parse_target(arg3)

# -----------------------------------------------------
# Filepaths
# -----------------------------------------------------
rootfile_dir = "/work/hallc/c-rsidis/skimfiles/pass0"
mc_dir = "/work/hallc/c-rsidis/relder/mc-single-arm"
input_settings_filepath = f"../../FILTER_type/{target_shortname.upper()}/{selected_run_type}_{selected_beam_pass}pass_{target_shortname}_runs.dat"
output_dir = f"{target_shortname.upper()}"

# -----------------------------------------------------
# Reading in monte-carlo report to obtain normfac
# -----------------------------------------------------
if target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
    mc_filepath = f"{mc_dir}/worksim/{selected_run_type}_{selected_beam_pass}pass_{target_shortname.lower()}.root"
    if not os.path.exists(mc_filepath):
        print(f"WARNING:\tNo mc-single-arm generated root file found: '{mc_filepath}'. Exiting...")
        exit(1)
    mc_report_filepath = f"{mc_dir}/outfiles/{selected_run_type}_{selected_beam_pass}pass_{target_shortname.lower()}.out"
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

with open(input_settings_filepath, "r", newline="") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        try:
            #runnums.append(int(row["runnum"]))
            runnums.append(row["runnum"])
            charge.append(row["qbeam"])
            weight.append(row["weight"])
            hms_p_val = float(row["hms_p"])
            polarity.append("-" if hms_p_val < 0 else "+")
        except KeyError:
            continue
        except ValueError:
            continue
print(f"Found {len(runnums)} runs, processing...")

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
    "H_dc_yp_fp": "hsypfp",
    "H_kin_W2": "w" #Altering the filling of this later
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
        "H_kin_W2": dict(binnum = 20, min = 9.6, max = 12.3),
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
        "H_kin_W2": dict(binnum = 20, min = 9.6, max = 12.3),
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

         if target_shortname == "dummy":

             hist_lh2 = bh.Histogram(axis, storage=bh.storage.Weight())
             hist_ld2 = bh.Histogram(axis, storage=bh.storage.Weight())

             yvals = df_cut["H_gtr_y"].values

             if var == "H_kin_W2":
                 if "H_kin_W" not in df_cut.columns:
                     raise KeyError("H_kin_W branch missing - cannot compute H_kin_W2 for dummy!")
                 varvals = df_cut["H_kin_W"].values**2
             else:
                 varvals = df_cut[var].values

             y_mid = 0.0
             upstream_mask = yvals < y_mid
             downstream_mask = yvals >= y_mid

             # Here is dummy:cell wall scaling if averaging across entire slab/cell wall
             # R_dummy_lh2_up = R_dummy_lh2_down = 1 / 7.2323
             # R_dummy_ld2_up = R_dummy_ld2_down = 1 / 7.7552

             # Here is dummy:cell wall scaling if considering separately upstream from downstream
             R_dummy_lh2_up, R_dummy_lh2_down = 1 / 8.2325, 1 / 6.4468
             R_dummy_ld2_up, R_dummy_ld2_down = 1 / 9.4990, 1 / 6.5493

             weights_lh2 = np.zeros(len(yvals))
             weights_lh2[upstream_mask] = R_dummy_lh2_up * run_weight
             weights_lh2[downstream_mask] = R_dummy_lh2_down * run_weight

             weights_ld2 = np.zeros(len(yvals))
             weights_ld2[upstream_mask] = R_dummy_ld2_up * run_weight
             weights_ld2[downstream_mask] = R_dummy_ld2_down * run_weight

             hist_lh2.fill(varvals, weight=weights_lh2)
             hist_ld2.fill(varvals, weight=weights_ld2)

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
             if var == "H_kin_W2":
                 values = df_cut["H_kin_W"].values**2
             else:
                 values = df_cut[var].values
                 
             hist.fill(values, weight=np.full(len(values), run_weight))

             counts = hist.view().value
             errors = np.sqrt(hist.view().variance)
         
             hist_data[var].append([runnum, charge[i], polarity[i]]+counts.tolist())
             hist_err_data[var].append([runnum, charge[i], polarity[i]]+errors.tolist())

# -----------------------------------------------------
# Monte carlo histogram and csv creation
# -----------------------------------------------------
if target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
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
    
    if target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        mc_var = variable_mc_map.get(var)
        
    if mc_var is None or target_shortname in {"dummy", "optics1", "optics2", "hole"}:
        continue
    
    axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"])
    if target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        hist_mc = bh.Histogram(axis, storage=bh.storage.Weight())

        if "weight" in df_mc_cut.columns:
            delta_temp = df_mc_cut["hsdelta"].values
            deltacorr = 1.0

            if target_shortname in {"al", "c", "cu"}:
                # Determined only with carbon,
                a = 1.012441e+00
                b = 3.055522e-03
                c = -1.111970e-03
                d = -6.311775e-05
                e = 1.411932e-05

                deltacorr = (a + b * delta_temp + c * delta_temp**2 + d * delta_temp**3 + e * delta_temp**4)

            elif target_shortname in {"ld2", "lh2", "dummy"}:
                a = 1.011192e+00
                b = 5.168480e-03
                c = -1.104189e-03
                d = -9.446273e-05
                e = 1.550629e-05

                deltacorr = (a + b * delta_temp + c * delta_temp**2 + d * delta_temp**3 + e * delta_temp**4)

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
        if var == "H_kin_W2":
            mc_values = df_mc_cut["w"].values**2
        else:
            mc_values = df_mc_cut[mc_var].values
            
        hist_mc.fill(mc_values, weight=event_weights)
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
    if target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        mc_row = ["MC", "0", "mc"] + mc_hist_data[var]
        hist_df = pd.concat([pd.DataFrame([mc_row], columns=columns), hist_df], ignore_index = True)
    
    output_filename = f"{selected_run_type}_{selected_beam_pass}pass_{target_shortname}_{var}_histo.csv"
    output_filepath = f"{output_dir}/{output_filename}"
    hist_df.to_csv(output_filepath, index=False)

    # Saving here for error bars
    hist_err_df = pd.DataFrame(hist_err_data[var], columns=columns)
    if target_shortname not in {"dummy", "optics1", "optics2", "hole"}:
        mc_err_row = ["MC", "0", "0"] + mc_hist_err[var]
        hist_err_df = pd.concat([pd.DataFrame([mc_err_row], columns=columns), hist_err_df], ignore_index=True)
    
    output_err_filename = f"{selected_run_type}_{selected_beam_pass}pass_{target_shortname}_{var}_err.csv"
    output_err_filepath = f"{output_dir}/{output_err_filename}"
    hist_err_df.to_csv(output_err_filepath, index=False)
    
print(f"Saved output CSV files to {target_shortname.upper()}/ folder.")
