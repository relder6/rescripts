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
from INIT.config import get_data_cuts, get_common_values
from INIT.config import parse_run_type, parse_beam_pass, parse_target, parse_bins

# -----------------------------------------------------
# Handling user inputs, listing directories
# -----------------------------------------------------
R_ld2 = 0.2557798343632005

arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None
arg3 = sys.argv[3] if len(sys.argv) > 3 else None
arg4 = sys.argv[4] if len(sys.argv) > 4 else None

selected_run_type = parse_run_type(arg1)
num_abbrev, num_longname, num_shortname, num_A, num_Z = parse_target(arg2)
denom_abbrev, denom_longname, denom_shortname, denom_A, denom_Z = parse_target(arg3)
nbins = parse_bins(arg4)

A_ratio = num_A / denom_A
num_N = num_A - num_Z
denom_N = denom_A - denom_Z

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
    csv_files.append(f"{xsec_ratio_dir}/XSEC_RATIO_{selected_run_type}_{beam_pass}_{num_abbrev}_to_{denom_abbrev}.csv")

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

    df["setting"] = f"RSIDIS({pass_label}Pass)"
    df["ebeam"] = ebeam
    df["theta"] = theta
    df["w2"] = df["w"] ** 2
    df["bc_corr"] = 0.0
    df["xsec_model_num"] = df["modelxsec_num"]
    df["xsec_model_denom"] = df["modelxsec_denom"]
    df["num_A"] = num_A
    df["num_Z"] = num_Z
    df["denom_A"] = denom_A
    df["denom_Z"] = denom_Z

    df_out = df[["setting", "num_A", "num_Z","denom_A","denom_Z","ebeam","theta","eprime","xbj", "q2", "w2", "epsilon", "xsec_exp_num", "xsec_exp_err_num", "xsec_exp_denom", "xsec_exp_err_denom", "xsec_ratio", "xsec_ratio_err", "bc_corr","xsec_model_num", "xsec_model_denom"]]

    df_out = df_out.dropna()

    all_rows.append(df_out)

df_data = pd.concat(all_rows, ignore_index=True)

output_csv = f"CSVs/DELTA_R_{selected_run_type.upper()}_bin_centered_{num_abbrev}_to_{denom_abbrev}.csv"

# -----------------------------------------------------
# Determining the bin centers
# -----------------------------------------------------
selected_var = "xbj"

df_xbj_minmax = (df_data.groupby("setting")[selected_var].agg(["min", "max"]).reset_index())
df_xbj_minmax = df_xbj_minmax.rename(columns={"min": f"{selected_var}_min", "max": f"{selected_var}_max"})
overlap_min = df_xbj_minmax[f"{selected_var}_min"].max()
overlap_max = df_xbj_minmax[f"{selected_var}_max"].min()
overlap_range = float(overlap_max) - float(overlap_min)

if nbins == 1:
    df_fit = df_data

    settings = df_fit["setting"].unique()
    if len(settings) < 2:
        raise ValueError(f"Need at least 2 experiments for common-center fit, got {len(settings)}")

    slopes = []
    intercepts = []
    for setting in settings:
        sub = df_fit[df_fit["setting"] == setting]
        x = sub[selected_var].to_numpy()
        y = sub["q2"].to_numpy()
        if len(x) < 2:
            raise ValueError(f"Not enough points to fit line for {setting}")
        m, b = np.polyfit(x, y, 1)
        slopes.append(m)
        intercepts.append(b)

    slopes = np.array(slopes)
    intercepts = np.array(intercepts)

    x_min = df_fit[selected_var].min()
    x_max = df_fit[selected_var].max()

    if len(settings) == 2:
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
    settings = df_data["setting"].unique()
    for setting in settings:
        sub = df_data[(df_data["setting"] == setting) & mask]
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

    infile_names = {"bc_num": f"{selected_run_type}_{beam_pass}_{num_abbrev}",
                    "bc_denom": f"{selected_run_type}_{beam_pass}_{denom_abbrev}"}

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
                
                model_eprime = np.nan
                model_theta = np.nan
                model_xbj = np.nan
                model_q2 = np.nan
                model_xsec = np.nan
            
                found_header = False

                for line in lines:
                    if "ub/GeV/sr" in line:
                        found_header = True
                        continue

                    if found_header:
                        parts = line.split()
                        if len(parts) >= 6:
                            try:
                                model_eprime = float(parts[0])
                                model_theta = float(parts[1])
                                model_xbj = float(parts[2])
                                model_q2 = float(parts[3])
                                model_xsec = float(parts[-1])
                                break
                            except:
                                pass
                        
                if np.isnan(model_xsec):
                    print(f"Failed to parse xsec for {beam_pass}, bin {bin_center}")

                results_tmp[f"{label}_xsec_model"] = model_xsec
                results_tmp[f"{label}_eprime"] = model_eprime
                results_tmp[f"{label}_theta"] = model_theta
                results_tmp[f"{label}_xbj"] = model_xbj
                results_tmp[f"{label}_q2"] = model_q2
                

            except Exception as e:
                print(f"Error running calc_dis_xsec: {e}")
                
        model_results.append({"ebeam": ebeam,
                              "bc_eprime": results_tmp.get("bc_num_eprime", np.nan),
                              "bc_theta": results_tmp.get("bc_num_theta", np.nan),
                              "bc_xbj": results_tmp.get("bc_num_xbj", np.nan),
                              "bc_q2": results_tmp.get("bc_num_q2", np.nan),
                              "bc_xsec_model_num": results_tmp.get("bc_num_xsec_model", np.nan),
                              "bc_xsec_model_denom": results_tmp.get("bc_denom_xsec_model", np.nan),})
        
df_centers = pd.DataFrame(model_results)
edges_centers = np.zeros(len(bin_centers)+1)

if len(bin_centers) == 1:
    # Only one bin, set edges around center
    delta = 0.5 * stepsize 
    edges_centers[0] = bin_centers[0] - delta
    edges_centers[1] = bin_centers[0] + delta
else:
    # Multiple bins: use midpoints as before
    edges_centers[1:-1] = (bin_centers[:-1] + bin_centers[1:]) / 2
    edges_centers[0]  = bin_centers[0] - (bin_centers[1] - bin_centers[0])/2
    edges_centers[-1] = bin_centers[-1] + (bin_centers[-1] - bin_centers[-2])/2

df_centers["bin_num"] = np.digitize(df_centers["bc_xbj"], edges_centers) - 1

print(f"\nXSEC Model at Bin Centers")
print(df_centers)

# -----------------------------------------------------
# Computing additional values for the dataframe
# -----------------------------------------------------
# df_data["bc_xsec_model_num"] = np.nan
# df_data["bc_xsec_model_denom"] = np.nan
# df_data["bc_xbj"] = np.nan
# df_data["bc_q2"] = np.nan
# df_data["bc_eprime"] = np.nan
# df_data["bc_theta"] = np.nan

df_centers_bin = (df_centers.groupby(["bin_num", "ebeam"], as_index=False)[["bc_xsec_model_num", "bc_xsec_model_denom", "bc_xbj", "bc_q2", "bc_eprime", "bc_theta"]].mean())

# for _, row in df_centers_bin.iterrows():
#     mask = (df_data["ebeam"] == row["ebeam"]) & (df_data["bin_num"] == row["bin_num"])
#     df_data.loc[mask, "bc_xsec_model_num"] = row["bc_xsec_model_num"]
#     df_data.loc[mask, "bc_xsec_model_denom"] = row["bc_xsec_model_denom"]
#     df_data.loc[mask, "bc_xbj"] = row["bc_xbj"]
#     df_data.loc[mask, "bc_q2"] = row["bc_q2"]
#     df_data.loc[mask, "bc_eprime"] = row["bc_eprime"]
#     df_data.loc[mask, "bc_theta"] = row["bc_theta"]

df_data = df_data.merge(df_centers_bin[["bin_num", "ebeam", "bc_xsec_model_num", "bc_xsec_model_denom", "bc_xbj", "bc_q2", "bc_eprime", "bc_theta"]], on=["bin_num", "ebeam"], how="left",)

df_data["bc_corr"] = ((df_data["bc_xsec_model_num"] / df_data["bc_xsec_model_denom"]) /(df_data["xsec_model_num"] / df_data["xsec_model_denom"]))

m_p = 0.93827208943 # proton mass, GeV

alpha = 1 / 137.035999177 # fine structure constant

df_data["theta_rad"] = np.deg2rad(df_data["theta"])

df_data["bc_nu"] = df_data["ebeam"] - df_data["bc_eprime"]

df_data["bc_epsilon"] = (1 + 2 * ( 1 + (df_data["bc_nu"]**2 / df_data["bc_q2"])) * np.tan(df_data["theta_rad"]/2)**2)**(-1)

df_data["bc_epsilon_p"] = df_data["bc_epsilon"] / ( 1 + df_data["bc_epsilon"] * R_ld2)

df_data["bc_gamma"] = alpha / (2 * np.pi**2 * df_data["bc_q2"]) * (df_data["ebeam"] - df_data["bc_nu"]) / (df_data["ebeam"]) * (df_data["bc_nu"] * (1 - df_data["bc_xbj"])) / (1 - df_data["bc_epsilon"])

df_data["sigma_num_to_sigma_denom"] = df_data["xsec_ratio"]

df_data["sigma_num_to_sigma_denom_err"] = df_data["xsec_ratio_err"]

df_data["bc_sigma_num_to_sigma_denom"] = df_data["xsec_ratio"] * df_data["bc_corr"]

df_data["bc_sigma_num_to_sigma_denom_err"] = df_data["xsec_ratio_err"] * df_data["bc_corr"]

df_data["epsilon_p"] = df_data["epsilon"] / (1 + df_data["epsilon"] * R_ld2)

col_final = ["setting", "num_A", "num_Z", "denom_A", "denom_Z", "ebeam", "theta", "theta_rad", "eprime", "xbj",
             "q2", "w2", "epsilon", "epsilon_p", "xsec_exp_num", "xsec_exp_err_num", "xsec_exp_denom", "xsec_exp_err_denom",
             "xsec_ratio", "xsec_ratio_err", "sigma_num_to_sigma_denom", "sigma_num_to_sigma_denom_err",
             "bin_num", "bc_xsec_model_num", "bc_xsec_model_denom", "bc_corr",
             "bc_xbj", "bc_q2", "bc_nu", "bc_epsilon", "bc_epsilon_p", "bc_gamma",
             "bc_sigma_num_to_sigma_denom", "bc_sigma_num_to_sigma_denom_err"]

df_final = df_data[col_final]

df_final.to_csv(output_csv, index=False)

print(f"Saved compiled CSV → {output_csv}")
print(f"bc_corr range: {df_final['bc_corr'].min():.6f} → {df_final['bc_corr'].max():.6f}")

# df_data["bc_xbj"] = bin_centers[df_data["bin_num"].to_numpy()]

# df_data["bc_q2"] = df_data["q2"] * df_data["bc_xbj"] / df_data["xbj"]

# df_data["emc"] = slac_emc_fit(df_data["xbj"],df_data["num_A"])

# df_data["bc_emc"] = slac_emc_fit(df_data["bc_xbj"],df_data["num_A"])

# df_data["bc_corr"] = df_data["bc_emc"] / df_data["emc"]

# m_p = 0.93827208943 # proton mass, GeV

# alpha = 1 / 137.035999177 # fine structure constant

# df_data["bc_nu"] = (1 / (2 * m_p) ) * df_data["bc_q2"] / df_data["bc_xbj"]

# df_data["theta_rad"] = np.deg2rad(df_data["theta"])

# df_data["bc_epsilon"] = (1 + 2 * ( 1 + (df_data["bc_nu"]**2 / df_data["bc_q2"])) * np.tan(df_data["theta_rad"]/2)**2)**(-1)

# R_ld2 = 0.1725

# df_data["bc_epsilon_p"] = df_data["bc_epsilon"] / ( 1 + df_data["bc_epsilon"] * R_ld2)

# df_data["bc_gamma"] = alpha / (2 * np.pi**2 * df_data["bc_q2"]) * (df_data["ebeam"] - df_data["bc_nu"]) / (df_data["ebeam"]) * (df_data["bc_nu"] * (1 - df_data["bc_xbj"])) / (1 - df_data["bc_epsilon"])

# df_data["bc_sigma_num_to_sigma_denom"] = df_data["xsec_ratio"] * df_data["bc_corr"]

# df_data["bc_sigma_num_to_sigma_denom_err"] = df_data["xsec_ratio_err"] * df_data["bc_corr"]

# col_final = ["setting", "num_A", "num_Z", "denom_A", "denom_Z", "ebeam", "theta", "theta_rad", "eprime", "xbj",
#              "q2", "w2", "epsilon", "xsec_exp_num", "xsec_exp_err_num", "xsec_exp_denom", "xsec_exp_err_denom",
#              "xsec_ratio", "xsec_ratio_err",
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
plt.title(f"{num_longname}/{denom_longname} Binning Test (nbins={nbins})")

plt.legend()
plt.grid(axis="both", linestyle="--", alpha=0.8)

plt.savefig(f"PNGs/{selected_run_type}_{num_abbrev}_to_{denom_abbrev}_binning_test.png")

#plt.show()

plt.close()
    
