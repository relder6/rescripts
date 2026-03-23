#!/usr/bin/env python3

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import boost_histogram as bh
import os, re, sys
import subprocess
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import get_common_run_inputs, get_data_cuts, get_common_values
from INIT.slac_emc_fit import slac_emc_fit

# -----------------------------------------------------
# Handling user inputs, listing directories
# -----------------------------------------------------
if len(sys.argv) == 6:
    selected_run_type = sys.argv[1].strip().lower()
    selected_beam_pass = sys.argv[2].strip()
    selected_num = sys.argv[3].strip().lower()
    selected_denom = sys.argv[4].strip().lower()
    nbins = int(sys.argv[5])
else:
    selected_run_type = input("Enter desired run type (default HMSDIS): ").strip().lower()
    if not selected_run_type:
        selected_run_type = "hmsdis"
    selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ").strip()
    selected_num = input("Enter numerator target: ").strip().lower()
    selected_denom = input("Enter denominator target: ").strip().lower()
    nbins = int(input("Enter bin number: "))
    
selected_beam_pass_to_energy_prefix = {"1": "2.", "2": "4.", "3": "6.", "4": "8.", "5": "10."}

beam_prefix = selected_beam_pass_to_energy_prefix.get(selected_beam_pass)

if not beam_prefix:
    print(f"Unknown pass: {selected_beam_pass}.  Please try again.")
    exit(1)

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

# Targets
num_short = selected_target_shortcut_to_target_variable.get(selected_num)
if not num_short:
    print(f"Unknown target: {selected_num}. Please try again.")
    exit(1)
num_long = selected_target_shortname_to_title_longname[num_short]
    
denom_short = selected_target_shortcut_to_target_variable.get(selected_denom)
if not denom_short:
    print(f"Unknown target: {selected_denom}. Please try again.")
    exit(1)
denom_long = selected_target_shortname_to_title_longname[denom_short]

selected_target_shortname_to_AZ = {"al":   (27, 13),
                                   "c":    (12, 6),
                                   "cu":   (64, 29),
                                   "ld2":  (2, 1),
                                   "lh2":  (1, 1)}

A_num, Z_num = selected_target_shortname_to_AZ.get(num_short)

N_num = A_num - Z_num

A_denom, Z_denom = selected_target_shortname_to_AZ.get(denom_short, (0.0, 0.0))

A_ratio = A_num / A_denom

vals = get_common_values()
ebeam_4pass = vals["ebeam_4pass"]
theta_4pass = vals["angle_4pass"]
ebeam_5pass = vals["ebeam_5pass"]
theta_5pass = vals["angle_5pass"]

model_xsec_dir = "../../../mc-single-arm/util/dis_xec"

xsec_ratio_dir = "../../XSEC/FORM_xsec/RATIOS"

beam_passes = ["4pass", "5pass"]

# -----------------------------------------------------
# Collecting extracted XSECs
# -----------------------------------------------------
csv_files = []

for beam_pass in beam_passes:
    csv_files.append(f"{xsec_ratio_dir}/XSEC_RATIO_{selected_run_type}_{beam_pass}_{num_short}_to_{denom_short}.csv")

all_rows = []

for filepath in csv_files:
    if not os.path.exists(filepath):
        print(f"Skipping missing file: {filepath}")
        continue

    df = pd.read_csv(filepath)

    m = re.search(r'_(\d)pass_', os.path.basename(filepath))
    if m:
        pass_label = m.group(1)
    else:
        print(f"Could not determine pass from filename: {filepath}, defaulting to unknown.")
        pass_label = "???"

    if int(pass_label) == 4:
        ebeam = ebeam_4pass
        theta = theta_4pass
    elif int(pass_label) == 5:
        ebeam = ebeam_5pass
        theta = theta_5pass

    # A,Z,eprime,theta,xbj,q2,w,epsilon,modelxsec,xsec_exp,xsec_exp_err

    df["exp"] = f"RSIDIS({pass_label}Pass)"
    df["ebeam"] = ebeam
    df["theta"] = theta
    df["w2"] = df["w"] ** 2
    df["bc_corr"] = 0.0
    df["xsec_model_num"] = df["modelxsec_num"]
    df["xsec_model_denom"] = df["modelxsec_denom"]

    df_out = df[["exp", "A_num", "Z_num","A_denom","Z_denom","ebeam","theta","eprime","xbj", "q2", "w2", "epsilon", "xsec_exp_num", "xsec_exp_err_num", "xsec_exp_denom", "xsec_exp_err_denom", "xsec_ratio_final", "xsec_ratio_final_err", "bc_corr","xsec_model_num", "xsec_model_denom"]]

    df_out = df_out.dropna()

    all_rows.append(df_out)

df_data = pd.concat(all_rows, ignore_index=True)

output_csv = f"CSVs/DELTA_R_{selected_run_type.upper()}_bin_centered_{num_short}_to_{denom_short}.csv"

# -----------------------------------------------------
# Determining the bin centers
# -----------------------------------------------------
selected_var = "xbj"

df_xbj_minmax = (df_data.groupby("exp")[selected_var].agg(["min", "max"]).reset_index())
df_xbj_minmax = df_xbj_minmax.rename(columns={"min": f"{selected_var}_min", "max": f"{selected_var}_max"})
overlap_min = df_xbj_minmax[f"{selected_var}_min"].max()
overlap_max = df_xbj_minmax[f"{selected_var}_max"].min()
overlap_range = float(overlap_max) - float(overlap_min)

if nbins == 1:
    df_fit = df_data

    exps = df_fit["exp"].unique()
    if len(exps) < 2:
        raise ValueError(f"Need at least 2 experiments for common-center fit, got {len(exps)}")

    slopes = []
    intercepts = []
    for exp in exps:
        sub = df_fit[df_fit["exp"] == exp]
        x = sub[selected_var].to_numpy()
        y = sub["q2"].to_numpy()
        if len(x) < 2:
            raise ValueError(f"Not enough points to fit line for {exp}")
        m, b = np.polyfit(x, y, 1)
        slopes.append(m)
        intercepts.append(b)

    slopes = np.array(slopes)
    intercepts = np.array(intercepts)

    x_min = df_fit[selected_var].min()
    x_max = df_fit[selected_var].max()

    if len(exps) == 2:
        slope1, slope2 = slopes
        int1, int2 = intercepts
        if np.isclose(slope1, slope2):
            x_center = 0.5 * (x_min + x_max)
        else:
            x_center = (int2 - int1) / (slope1 - slope2)
    else:
        x_grid = np.linspace(x_min, x_max, 1000)
        y_all = slopes[:, None] * x_grid[None, :] + intercepts[:, None]
        var_y = np.var(y_all, axis=0)
        idx_min = np.argmin(var_y)
        x_center = x_grid[idx_min]

    x_center = max(min(x_center, x_max), x_min)

    bin_centers = np.array([x_center])

    stepsize = x_max - x_min
    edges = np.array([x_min, x_max])

    df_data["bin_num"] = 0

if nbins != 1:
    stepsize = float(overlap_range) / float(nbins)
    central_value = overlap_min + overlap_range / 2

    if nbins % 2 == 1:
        central_value = (overlap_min + overlap_max) / 2
        half_n = nbins // 2
        bin_centers = central_value + stepsize * np.arange(-half_n, half_n + 1)
    else:
        central_value = (overlap_min + overlap_max) / 2
        half_n = nbins // 2
        bin_centers = central_value + stepsize * (np.arange(-half_n, half_n) + 0.5)

    bin_centers = np.array(sorted(bin_centers))

    edges = np.zeros(len(bin_centers) + 1)
    edges[1:-1] = (bin_centers[:-1] + bin_centers[1:]) / 2.0
    edges[0]  = overlap_min  - (bin_centers[1]  - bin_centers[0]) / 2.0
    edges[-1] = overlap_max + (bin_centers[-1] - bin_centers[-2]) / 2.0

    centers_arr = bin_centers.reshape(1, -1)
    vals = df_data[selected_var].to_numpy().reshape(-1, 1)
    df_data["bin_num"] = np.argmin(np.abs(vals - centers_arr), axis=1)

print(f"Overlap range: {overlap_min:.4f} → {overlap_max:.4f}")
print(f"Stepsize: {stepsize:.4f}")
print(f"Bin centers ({nbins} bins):\n{bin_centers}")

bin_info = []

for bin_idx, bin_center in enumerate(bin_centers):
    mask = df_data["bin_num"] == bin_idx
    bc_q2_list = []
    exps = df_data["exp"].unique()
    for exp in exps:
        sub = df_data[(df_data["exp"] == exp) & mask]
        if len(sub) < 2:
            continue
        x = sub[selected_var].to_numpy()
        y = sub["q2"].to_numpy()
        m, b = np.polyfit(x, y, 1)
        bc_q2 = m * bin_center + b
        bc_q2_list.append(bc_q2)
    avg_q2 = np.mean(bc_q2_list)
    bin_info.append({"bin_num": bin_idx, "bc_xbj": bin_center, "bc_q2": avg_q2})

df_bins = pd.DataFrame(bin_info)
print("Bin-centered x and Q² values:")
print(df_bins)

# -----------------------------------------------------
# Now building input strings, collecting model xsec of bin centers
# -----------------------------------------------------
model_results = []

for beam_pass in beam_passes:
    if beam_pass == "4pass":
        theta_inp = theta_4pass
        ebeam = ebeam_4pass
    if beam_pass == "5pass":
        theta_inp = theta_5pass
        ebeam = ebeam_5pass

    infile_names = {"num": f"{selected_run_type}_{beam_pass}_{num_short}",
                    "denom": f"{selected_run_type}_{beam_pass}_{denom_short}"}

    # The input string depends on the version of Dave's xsec tool, mc-single-arm/util/dis_xec/calc_dis_xsec
    # Right now, the input string is flag (0 = fixed theta, bin in eprime; 1 = fixed theta, bin in xbj; 2 = fixed Q2, bin in xbj...
    # fixed var (theta or Q2), <var>min, <var>step, <var>bin_num, then input filename

    for _, bin_row in df_bins.iterrows():
        bc_xbj = bin_row["bc_xbj"]
        bc_q2 = bin_row["bc_q2"]
        bin_num = bin_row["bin_num"]

        results_tmp = {}

        for label, infile_name in infile_names.items():

            input_string = f"2\n{bc_q2:.6f}\n{bc_xbj:.6f}\n1\n1\n{infile_name}\n"

            # print(f"INPUT STRING: {input_string}")

            try:
                result = subprocess.run(["./calc_dis_xsec"], input=input_string, text=True, capture_output=True, cwd=model_xsec_dir)
                # print("STDOUT:")
                # print(result.stdout)
                # print("STDERR:")
                # print(result.stderr)
                output = result.stdout

                lines = output.split("\n")

                model_xsec_value = np.nan
                model_q2_value = np.nan
            
                found_header = False

                for line in lines:
                    if "ub/GeV/sr" in line:
                        found_header = True
                        continue

                    if found_header:
                        parts = line.split()
                        if len(parts) >= 6:
                            try:
                                model_q2_value = float(parts[3])
                                model_xsec_value = float(parts[-1])
                                break
                            except:
                                pass
                        
                if np.isnan(model_xsec_value):
                    print(f"Failed to parse xsec for {beam_pass}, bin {bin_center}")

                results_tmp[label] = model_xsec_value
                results_tmp[f"{label}_q2"] = model_q2_value

            except Exception as e:
                print(f"Error running calc_dis_xsec: {e}")
                
        model_results.append({"ebeam": ebeam,
                              "bc_theta": theta_inp,
                              "bc_xbj": bin_center,
                              "bc_q2": results_tmp.get("num_q2", np.nan),
                              "xsec_model_num": results_tmp.get("num", np.nan),
                              "xsec_model_denom": results_tmp.get("denom", np.nan),})
        
df_centers = pd.DataFrame(model_results)
center_to_bin = dict(zip(bin_centers, range(len(bin_centers))))
df_centers["bin_num"] = df_centers["bc_xbj"].map(center_to_bin)

print(f"\nXSEC Model at Bin Centers")
print(df_centers)

df_centers_bin = (df_centers.groupby(["bin_num", "ebeam"], as_index=False)[["xsec_model_num", "xsec_model_denom", "bc_xbj", "bc_q2"]].mean())

df_data["bc_xsec_model_num"] = np.nan
df_data["bc_xsec_model_denom"] = np.nan
df_data["bc_xbj"] = np.nan
df_data["bc_q2"] = np.nan

for _, row in df_centers_bin.iterrows():
    mask = (df_data["ebeam"] == row["ebeam"]) & (df_data["bin_num"] == row["bin_num"])
    df_data.loc[mask, "bc_xsec_model_num"] = row["xsec_model_num"]
    df_data.loc[mask, "bc_xsec_model_denom"] = row["xsec_model_denom"]
    df_data.loc[mask, "bc_xbj"] = row["bc_xbj"]
    df_data.loc[mask, "bc_q2"] = row["bc_q2"]

df_data["bc_corr"] = ((df_data["bc_xsec_model_num"] / df_data["bc_xsec_model_denom"]) /(df_data["xsec_model_num"] / df_data["xsec_model_denom"]))

m_p = 0.93827208943 # proton mass, GeV

alpha = 1 / 137.035999177 # fine structure constant

R_ld2 = 0.2562

df_data["theta_rad"] = np.deg2rad(df_data["theta"])

df_data["bc_nu"] = (1 / (2 * m_p) ) * df_data["bc_q2"] / df_data["bc_xbj"]

df_data["bc_epsilon"] = (1 + 2 * ( 1 + (df_data["bc_nu"]**2 / df_data["bc_q2"])) * np.tan(df_data["theta_rad"]/2)**2)**(-1)

df_data["bc_epsilon_p"] = df_data["bc_epsilon"] / ( 1 + df_data["bc_epsilon"] * R_ld2)

df_data["bc_gamma"] = alpha / (2 * np.pi**2 * df_data["bc_q2"]) * (df_data["ebeam"] - df_data["bc_nu"]) / (df_data["ebeam"]) * (df_data["bc_nu"] * (1 - df_data["bc_xbj"])) / (1 - df_data["bc_epsilon"])

df_data["sigma_num_to_sigma_denom"] = df_data["xsec_ratio_final"]

df_data["sigma_num_to_sigma_denom_err"] = df_data["xsec_ratio_final_err"]

df_data["bc_sigma_num_to_sigma_denom"] = df_data["xsec_ratio_final"] * df_data["bc_corr"]

df_data["bc_sigma_num_to_sigma_denom_err"] = df_data["xsec_ratio_final_err"] * df_data["bc_corr"]

df_data["epsilon_p"] = df_data["epsilon"] / (1 + df_data["epsilon"] * R_ld2)

col_final = ["exp", "A_num", "Z_num", "A_denom", "Z_denom", "ebeam", "theta", "theta_rad", "eprime", "xbj",
             "q2", "w2", "epsilon", "epsilon_p", "xsec_exp_num", "xsec_exp_err_num", "xsec_exp_denom", "xsec_exp_err_denom",
             "xsec_ratio_final", "xsec_ratio_final_err", "sigma_num_to_sigma_denom", "sigma_num_to_sigma_denom_err",
             "bin_num", "bc_xsec_model_num", "bc_xsec_model_denom", "bc_corr",
             "bc_xbj", "bc_q2", "bc_nu", "bc_epsilon", "bc_epsilon_p", "bc_gamma",
             "bc_sigma_num_to_sigma_denom", "bc_sigma_num_to_sigma_denom_err"]

df_final = df_data[col_final]

df_final.to_csv(output_csv, index=False)

print(f"Saved compiled CSV → {output_csv}")
print(f"bc_corr range: {df_final['bc_corr'].min():.6f} → {df_final['bc_corr'].max():.6f}")

# df_data["bc_xbj"] = bin_centers[df_data["bin_num"].to_numpy()]

# df_data["bc_q2"] = df_data["q2"] * df_data["bc_xbj"] / df_data["xbj"]

# df_data["emc"] = slac_emc_fit(df_data["xbj"],df_data["A_num"])

# df_data["bc_emc"] = slac_emc_fit(df_data["bc_xbj"],df_data["A_num"])

# df_data["bc_corr"] = df_data["bc_emc"] / df_data["emc"]

# m_p = 0.93827208943 # proton mass, GeV

# alpha = 1 / 137.035999177 # fine structure constant

# df_data["bc_nu"] = (1 / (2 * m_p) ) * df_data["bc_q2"] / df_data["bc_xbj"]

# df_data["theta_rad"] = np.deg2rad(df_data["theta"])

# df_data["bc_epsilon"] = (1 + 2 * ( 1 + (df_data["bc_nu"]**2 / df_data["bc_q2"])) * np.tan(df_data["theta_rad"]/2)**2)**(-1)

# R_ld2 = 0.1725

# df_data["bc_epsilon_p"] = df_data["bc_epsilon"] / ( 1 + df_data["bc_epsilon"] * R_ld2)

# df_data["bc_gamma"] = alpha / (2 * np.pi**2 * df_data["bc_q2"]) * (df_data["ebeam"] - df_data["bc_nu"]) / (df_data["ebeam"]) * (df_data["bc_nu"] * (1 - df_data["bc_xbj"])) / (1 - df_data["bc_epsilon"])

# df_data["bc_sigma_num_to_sigma_denom"] = df_data["xsec_ratio_final"] * df_data["bc_corr"]

# df_data["bc_sigma_num_to_sigma_denom_err"] = df_data["xsec_ratio_final_err"] * df_data["bc_corr"]

# col_final = ["exp", "A_num", "Z_num", "A_denom", "Z_denom", "ebeam", "theta", "theta_rad", "eprime", "xbj",
#              "q2", "w2", "epsilon", "xsec_exp_num", "xsec_exp_err_num", "xsec_exp_denom", "xsec_exp_err_denom",
#              "xsec_ratio_final", "xsec_ratio_final_err",
#              "emc", "bin_num", "bc_emc", "bc_corr",
#              "bc_xbj", "bc_q2", "bc_nu", "bc_epsilon", "bc_epsilon_p", "bc_gamma",
#              "bc_sigma_num_to_sigma_denom", "bc_sigma_num_to_sigma_denom_err"]

# df_final = df_data[col_final]

# df_final.to_csv(output_csv, index=False)

# print(f"Saved compiled CSV → {output_csv}")
# print(f"bc_corr range: {df_final['bc_corr'].min():.6f} → {df_final['bc_corr'].max():.6f}")

# -----------------------------------------------------
# Plotting
# -----------------------------------------------------
plt.figure()

unique_bins = sorted(df_data["bin_num"].unique())

if nbins != 1:
    plt.axvspan(overlap_min, overlap_max, color="lightskyblue", alpha=0.1, label="Overlap Region")
else:
    plt.axvspan(x_min, x_max, color="lightskyblue", alpha=0.1, label="Data Range")

plt.scatter(df_data["bc_xbj"],df_data["bc_q2"], marker="*",label="Bin centers", s=64)

for bin_num in unique_bins:
    mask = df_data["bin_num"] == bin_num
    plt.scatter(df_data.loc[mask, "xbj"],df_data.loc[mask, "q2"], label=f"Bin {bin_num}", s=12)

for i in range(1, len(edges)-1):
    if i == 1:
        plt.axvline(edges[i], ls="--", color="red", alpha=0.8, label="Bin Edges")
    else:
        plt.axvline(edges[i], ls="--", color="red", alpha=0.8)

plt.xlabel(r"x$_{bj}$")
plt.ylabel(r"Q$^2$")
plt.title(f"{num_long}/{denom_long} Binning Test (nbins={nbins})")

plt.legend()
plt.grid(axis="both", linestyle="--", alpha=0.8)

plt.savefig(f"PNGs/{selected_run_type}_{num_short}_to_{denom_short}_binning_test.png")

plt.show()

plt.close()
    
