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
from INIT.config import parse_run_type, parse_beam_pass, parse_target, get_data_cuts, get_flags, parse_phase

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
rootfile_type = 1 #0 for full rootfile, 1 for skimfiles

flags = get_flags()

USING_DELTA_CORR = flags["USING_DELTA_CORR"]

arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None
arg3 = sys.argv[3] if len(sys.argv) > 3 else None
arg4 = sys.argv[4] if len(sys.argv) > 4 else None

selected_run_type = parse_run_type(arg1)
selected_beam_pass, beam_prefix = parse_beam_pass(arg2)
target_abbrev, target_longname, target_shortname, target_A, target_Z = parse_target(arg3)
phase = parse_phase(arg4)

# -----------------------------------------------------
# Filepaths
# -----------------------------------------------------
if phase == "I":
    rootfile_dir = "/work/hallc/c-rsidis/skimfiles/pass0p1"
elif phase == "II":
    rootfile_dir = "/volatile/hallc/c-rsidis/relder/STUFF"
mc_dir = "/work/hallc/c-rsidis/relder/mc-single-arm"
if target_abbrev not in {"dummy_up", "dummy_down"}:
    input_settings_filepath = f"../../FILTER_type/{target_abbrev.upper()}/filtered_{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}.csv"
else:
    input_settings_filepath = f"../../FILTER_type/DUMMY/filtered_{selected_run_type}_{selected_beam_pass}pass_phase{phase}_dummy.csv"
output_dir = f"{target_abbrev.upper()}"

if not os.path.exists(input_settings_filepath):
    print(f"File {input_settings_filepath} not found; Exiting...")
    sys.exit(0)

# -----------------------------------------------------
# Reading in monte-carlo report to obtain normfac
# -----------------------------------------------------
mc_exists = False
normfac = 0.0

if target_abbrev not in {"dummy", "optics1", "optics2", "hole"}:
    mc_filepath = f"{mc_dir}/worksim/{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev.lower()}.root"

    mc_exists = os.path.exists(mc_filepath)

    if mc_exists:
        mc_report_filepath = f"{mc_dir}/outfiles/{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev.lower()}.out"
        with open(mc_report_filepath, "r") as infile:
            normfac_line = [line for line in infile if "NORMFAC" in line.upper()]
        normfac = float(normfac_line[0].split(":")[1].split()[0]) if normfac_line else None
    else:
        print(f"WARNING:\tNo mc-single-arm generated root file found: '{mc_filepath}'. Zeroes will be written in MC row.")    
    # print(f"DEBUG: normfac found:{normfac}")

# -----------------------------------------------------
# Defining branches, using uproot to put them in data frames
# -----------------------------------------------------
runnums = []
weight = []
charge = []
current = []
polarity = []

with open(input_settings_filepath, "r", newline="") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        try:
            #runnums.append(int(row["runnum"]))
            runnums.append(row["runnum"])
            charge.append(row["qbeam_2"])
            current.append(row["ibeam_2"])
            weight.append(row["weight"])
            hms_p_val = float(row["hms_p"])
            polarity.append("-" if hms_p_val < 0 else "+")
        except KeyError:
            continue
        except ValueError:
            continue

# -----------------------------------------------------
# Defining branches, using uproot to put them in data frames
# -----------------------------------------------------
branches_mc = ["hsdelta", "q2", "xb", "w", "weight", "eprime", "hsytar", "hsxptar", "hsyptar", "hsxfp", "hsxpfp", "hsyfp", "hsypfp"]

branches = ["H.gtr.dp", "H.cal.etottracknorm", "H.gtr.ph",
            "H.gtr.th", "H.gtr.x", "H.gtr.y",
            "H.kin.Q2", "H.kin.x_bj", "H.kin.W",
            "H.cer.npeSum", "H.gtr.p",
            "H.dc.x_fp", "H.dc.xp_fp", "H.dc.y_fp", "H.dc.yp_fp"]

variable_mc_map = {"H.gtr.dp": "hsdelta",
                   "H.gtr.ph": "hsyptar",
                   "H.gtr.th": "hsxptar",
                   "H.kin.Q2": "q2",
                   "H.kin.x_bj": "xb",
                   "H.kin.W": "w",
                   "H.gtr.p": "eprime",
                   "H.gtr.y": "hsytar",
                   "H.gtr.th": "hsxptar",
                   "H.gtr.ph": "hsyptar",
                   "H.dc.x_fp": "hsxfp",
                   "H.dc.xp_fp": "hsxpfp",
                   "H.dc.y_fp": "hsyfp",
                   "H.dc.yp_fp": "hsypfp",
                   "H.kin.W2": "w"}

if rootfile_type == 1:
    branches = [branch.replace(".", "_") for branch in branches]
    variable_mc_map = {key.replace(".", "_"): value for key, value in variable_mc_map.items()}

# -----------------------------------------------------
# Binning
# -----------------------------------------------------
if selected_beam_pass == "3":
    custom_bins = {"H.gtr.dp": dict(binnum = 20, min = -10.000, max = 10.000),
                   "H.gtr.ph": dict(binnum = 20, min = -0.050, max = 0.050),
                   "H.gtr.th": dict(binnum = 20, min = -0.100, max = 0.100),
                   "H.kin.Q2": dict(binnum = 20, min = 2.900, max = 6.000),
                   "H.kin.x_bj": dict(binnum = 20, min = 0.2, max = 0.7),
                   "H.kin.W": dict(binnum = 20, min = 2.000, max = 3.000),
                   "H.gtr.p": dict(binnum = 100, min = 1.000, max = 1.400),
                   "H.gtr.y": dict(binnum = 100, min = -4.0, max = 4.0),
                   "H.gtr.th": dict(binnum = 100, min = -0.1, max = 0.1),
                   "H.gtr.ph": dict(binnum = 100, min = -0.05, max = 0.05),
                   "H.dc.x_fp": dict(binnum = 20, min = -50, max = 50),
                   "H.dc.xp_fp": dict(binnum = 20, min = -0.08, max = 0.08),
                   "H.dc.y_fp": dict(binnum = 20, min = -30, max = 30),
                   "H.dc.yp_fp": dict(binnum = 20, min = -0.04, max = 0.04),
                   "H.kin.W2": dict(binnum = 20, min = 4.00, max = 9.00),}
    
if selected_beam_pass == "4":
    custom_bins = {"H.gtr.dp": dict(binnum = 20, min = -10.000, max = 10.000),
                   "H.gtr.ph": dict(binnum = 20, min = -0.050, max = 0.050),
                   "H.gtr.th": dict(binnum = 20, min = -0.100, max = 0.100),
                   "H.kin.Q2": dict(binnum = 20, min = 2.400, max = 4.200),
                   "H.kin.x_bj": dict(binnum = 20, min = 0.175, max = 0.325),
                   "H.kin.W": dict(binnum = 20, min = 3.100, max = 3.500),
                   "H.gtr.p": dict(binnum = 100, min = 1.3, max = 1.8),
                   "H.gtr.y": dict(binnum = 100, min = -4.0, max = 4.0),
                   "H.gtr.th": dict(binnum = 100, min = -0.1, max = 0.1),
                   "H.gtr.ph": dict(binnum = 100, min = -0.05, max = 0.05),
                   "H.dc.x_fp": dict(binnum = 20, min = -50, max = 50),
                   "H.dc.xp_fp": dict(binnum = 20, min = -0.08, max = 0.08),
                   "H.dc.y_fp": dict(binnum = 20, min = -30, max = 30),
                   "H.dc.yp_fp": dict(binnum = 20, min = -0.04, max = 0.04),
                   "H.kin.W2": dict(binnum = 20, min = 9.6, max = 12.3),}
    
if selected_beam_pass == "5":
    custom_bins = {"H.gtr.dp": dict(binnum = 20, min = -10.000, max = 10.000),
                   "H.gtr.ph": dict(binnum = 20, min = -0.050, max = 0.050),
                   "H.gtr.th": dict(binnum = 20, min = -0.100, max = 0.100),
                   "H.kin.Q2": dict(binnum = 20, min = 2.400, max = 4.200),
                   "H.kin.x_bj": dict(binnum = 20, min = 0.175, max = 0.325),
                   "H.kin.W": dict(binnum = 20, min = 3.100, max = 3.500),
                   "H.gtr.p": dict(binnum = 100, min = 3.2, max = 5.2),
                   "H.gtr.y": dict(binnum = 100, min = -4.0, max = 4.0),
                   "H.gtr.th": dict(binnum = 100, min = -0.1, max = 0.1),
                   "H.gtr.ph": dict(binnum = 100, min = -0.05, max = 0.05),
                   "H.dc.x_fp": dict(binnum = 20, min = -50, max = 50),
                   "H.dc.xp_fp": dict(binnum = 20, min = -0.08, max = 0.08),
                   "H.dc.y_fp": dict(binnum = 20, min = -30, max = 30),
                   "H.dc.yp_fp": dict(binnum = 20, min = -0.04, max = 0.04),
                   "H.kin.W2": dict(binnum = 20, min = 9.6, max = 12.3),}

if rootfile_type == 1:
    custom_bins = {key.replace(".", "_"): value for key, value in custom_bins.items()}
    

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
     if rootfile_type == 0:
         rootfile_path = f"{rootfile_dir}/hms_coin_replay_production_{runnum}_-1.root"
     if rootfile_type == 1:
         rootfile_path = f"{rootfile_dir}/skimmed_hms_coin_replay_production_{runnum}_-1.root"
     if not os.path.exists(rootfile_path):
         print(f"WARNING: Missing {rootfile_path}, skipping...")
         continue
     tree = uproot.open(rootfile_path)["T"]
     arr = tree.arrays(branches, library = "np")
     df = pd.DataFrame(arr)
     
     if rootfile_type == 0:
         dp = df["H.gtr.dp"]
         cer = df["H.cer.npeSum"]
         cal = df["H.cal.etottracknorm"]
         ytar = df["H.gtr.y"]
     elif rootfile_type == 1:
         dp = df["H_gtr_dp"]
         cer = df["H_cer_npeSum"]
         cal = df["H_cal_etottracknorm"]
         ytar = df["H_gtr_y"]
         
     if target_abbrev not in {"dummy_up", "dummy_down"}:

         data_cut = (dp.between(cuts["H_gtr_dp_min_cut"], cuts["H_gtr_dp_max_cut"]) &
                     (cer > cuts["H_cer_npeSum_cut"]) &
                     (cal > cuts["H_cal_etottracknorm_cut"]))
     elif target_abbrev == "dummy_up":
         data_cut = (dp.between(cuts["H_gtr_dp_min_cut"], cuts["H_gtr_dp_max_cut"]) &
                     (cer > cuts["H_cer_npeSum_cut"]) &
                     (cal > cuts["H_cal_etottracknorm_cut"]) &
                     (ytar < 0))
     elif target_abbrev == "dummy_down":
         data_cut = (dp.between(cuts["H_gtr_dp_min_cut"], cuts["H_gtr_dp_max_cut"]) &
                     (cer > cuts["H_cer_npeSum_cut"]) &
                     (cal > cuts["H_cal_etottracknorm_cut"]) &
                     (ytar >= 0))
     df_cut = df[data_cut].copy()

     run_weight = float(weight[i])

     for var, bins in custom_bins.items():
         axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"], underflow=True, overflow=True)

         if target_abbrev == "dummy":
             hist_lh2 = bh.Histogram(axis, storage=bh.storage.Weight())
             hist_ld2 = bh.Histogram(axis, storage=bh.storage.Weight())

             if rootfile_type == 0:
                 yvals = df_cut["H.gtr.y"].values
                 if var == "H.kin.W2":
                     if "H.kin.W" not in df_cut.columns:
                         raise KeyError("H.kin.W branch missing - cannot compute H.kin.W2 for dummy!")
                     varvals = df_cut["H.kin.W"].values**2
                 else:
                     varvals = df_cut[var].values
                     
             elif rootfile_type == 1:
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

             if phase == "I":
                 R_dummy_lh2_up, R_dummy_lh2_down = 1 / 8.2325, 1 / 6.4468
                 R_dummy_ld2_up, R_dummy_ld2_down = 1 / 9.4990, 1 / 6.5493
             elif phase == "II":
                 R_dummy_lh2_up, R_dummy_lh2_down = 1 / 4.43654, 1 / 3.48419
                 R_dummy_ld2_up, R_dummy_ld2_down = 1 / 5.11908, 1 / 3.53979

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

             hist_data.setdefault(var + "_lh2", []).append([runnum, charge[i], current[i], polarity[i]] + counts_lh2.tolist())
             hist_err_data.setdefault(var + "_lh2", []).append([runnum, charge[i], current[i], polarity[i]] + errors_lh2.tolist())

             hist_data.setdefault(var + "_ld2", []).append([runnum, charge[i], current[i], polarity[i]] + counts_ld2.tolist())
             hist_err_data.setdefault(var + "_ld2", []).append([runnum, charge[i], current[i], polarity[i]] + errors_ld2.tolist())
             
         else:     
             hist = bh.Histogram(axis, storage=bh.storage.Weight())
             if rootfile_type == 0:
                 if var == "H.kin.W2":
                     values = df_cut["H.kin.W"].values**2
                 else:
                     values = df_cut[var].values
             elif rootfile_type == 1:
                 if var == "H_kin_W2":
                     values = df_cut["H_kin_W"].values**2
                 else:
                     values = df_cut[var].values
                 
             hist.fill(values, weight=np.full(len(values), run_weight))

             counts = hist.view().value
             errors = np.sqrt(hist.view().variance)
         
             hist_data[var].append([runnum, charge[i], current[i], polarity[i]]+counts.tolist())
             hist_err_data[var].append([runnum, charge[i], current[i], polarity[i]]+errors.tolist())

# -----------------------------------------------------
# Monte carlo histogram and csv creation
# -----------------------------------------------------
mc_hist_data = {}
mc_hist_err = {}

if mc_exists:
    if target_abbrev not in {"dummy", "optics1", "optics2", "hole"}:

        mc_file = uproot.open(mc_filepath)
        mc_tree = mc_file["h10"]
        df_mc = pd.DataFrame(mc_tree.arrays(branches_mc, library="np"))
        mc_cut = (df_mc["hsdelta"].between(cuts["H_gtr_dp_min_cut"], cuts["H_gtr_dp_max_cut"]))
        df_mc_cut = df_mc[mc_cut].copy()
        if df_mc_cut.empty:
            print("ERROR: df_mc_cut is EMPTY — hsdelta cut removed everything.")
        mc_hist_data = {}
        mc_hist_err = {}

    for var, bins in custom_bins.items():
        mc_var = None

        if target_abbrev not in {"optics1", "optics2", "hole"}:
            mc_var = variable_mc_map.get(var)

        if mc_var is None or target_abbrev in {"dummy", "optics1", "optics2", "hole"}:
            continue

        axis = bh.axis.Regular(bins["binnum"], bins["min"], bins["max"])
        if target_abbrev not in {"dummy", "optics1", "optics2", "hole"}:
            hist_mc = bh.Histogram(axis, storage=bh.storage.Weight())

            if "weight" in df_mc_cut.columns:
                delta_temp = df_mc_cut["hsdelta"].values
                deltacorr = 1.0

                if target_abbrev in {"al", "c", "cu", "dummy_up", "dummy_down"}:
                    if phase == "I":
                        # Determined only with carbon, OLD
                        a = 1.012441e+00
                        b = 3.055522e-03
                        c = -1.111970e-03
                        d = -6.311775e-05
                        e = 1.411932e-05

                    elif phase == "II":
                        a = 1
                        b = 0
                        c = 0
                        d = 0
                        e = 0

                    deltacorr = (a + b * delta_temp + c * delta_temp**2 + d * delta_temp**3 + e * delta_temp**4)
                elif target_abbrev in {"ld2", "lh2", "dummy"}:
                    if phase == "I":
                        # Determined with both ld2, lh2, OLD
                        a = 1.011192e+00
                        b = 5.168480e-03
                        c = -1.104189e-03
                        d = -9.446273e-05
                        e = 1.550629e-05

                    elif phase == "II":
                        a = 1
                        b = 0
                        c = 0
                        d = 0
                        e = 0

                    deltacorr = (a + b * delta_temp + c * delta_temp**2 + d * delta_temp**3 + e * delta_temp**4)

                if USING_DELTA_CORR:
                    event_weights = df_mc_cut["weight"].values * normfac * deltacorr
                else:
                    event_weights = df_mc_cut["weight"].values * normfac
            else:
                print("Weights branch not found.  Exiting...")
                print(df_mc_cut.columns.tolist())
                exit(1)
            if var == "H_kin_W2":
                mc_values = df_mc_cut["w"].values**2
            else:
                mc_values = df_mc_cut[mc_var].values

            hist_mc.fill(mc_values, weight=event_weights)

            mc_hist_data[var] = hist_mc.view().value.tolist()
            mc_hist_err[var] = np.sqrt(hist_mc.view().variance).tolist()


# -----------------------------------------------------
# Saving consolidated CSV
# -----------------------------------------------------

all_rows = []

for var, rows in hist_data.items():

    base_var = var.replace("_lh2", "").replace("_ld2", "")
    output_var = base_var.replace(".", "_")

    bin_edges = bin_edges_dict[base_var]
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    nbins = len(bin_centers)
    bin_min = bin_edges[0]
    bin_max = bin_edges[-1]

    if target_abbrev not in {"dummy", "optics1", "optics2", "hole"}:

        if var in mc_hist_data:
            mc_values = mc_hist_data[var]
        else:
            mc_values = [0.0] * nbins

        all_rows.append([output_var, "mc", 0, 0, 0, "-",nbins, bin_min, bin_max] + mc_values)

        if var in mc_hist_err:
            mc_errors = mc_hist_err[var]
        else:
            mc_errors = [0.0] * nbins

        all_rows.append([output_var, "mc_err", 0, 0, 0, "-",nbins, bin_min, bin_max] + mc_errors)

    for row in rows:

        runnum = row[0]
        charge_val = row[1]
        current_val = row[2]
        polarity_val = row[3]
        counts = row[4:]

        all_rows.append([output_var, "data",runnum,charge_val,current_val,polarity_val,nbins,bin_min,bin_max] + counts)

    for row in hist_err_data[var]:

        runnum = row[0]
        charge_val = row[1]
        current_val = row[2]
        polarity_val = row[3]
        errors = row[4:]

        all_rows.append([output_var, "err",runnum,charge_val,current_val,polarity_val,nbins,bin_min,bin_max] + errors)

# Create column names
max_bins = max(len(row) - 9 for row in all_rows)

columns = ["variable","type","runnum","charge","current","polarity","nbins","bin_min","bin_max"]

columns += [f"bin{i}" for i in range(max_bins)]

output_df = pd.DataFrame(all_rows, columns=columns)

os.makedirs(output_dir, exist_ok=True)

output_filepath = (f"{output_dir}/{selected_run_type}_{selected_beam_pass}pass_phase{phase}_{target_abbrev}.csv")

output_df.to_csv(output_filepath, index=False)

print(f"Saved {len(runnums)} runs to CSV: {output_filepath}")
