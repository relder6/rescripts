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
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import get_common_run_inputs, get_data_cuts, get_common_values

# -----------------------------------------------------
# Handling user inputs, listing directories
# -----------------------------------------------------
selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, selected_target_titlename, selected_target_A, selected_target_Z = get_common_run_inputs()

vals = get_common_values()
ebeam_4pass = vals["ebeam_4pass"]
theta_4pass = vals["angle_4pass"]
ebeam_5pass = vals["ebeam_5pass"]
theta_5pass = vals["angle_5pass"]

nbins = int(input("Indicate number of bins: "))

model_xsec_dir = "../../mc-single-arm/util/dis_xec"

ratio_directory = "../XSEC/FORM_xsec"

beam_passes = ["4pass", "5pass"]

# -----------------------------------------------------
# Collecting extracted XSECs
# -----------------------------------------------------
csv_files = []

for beam_pass in beam_passes:
    csv_files.append(f"{ratio_directory}/{selected_target_shortname.upper()}/XSEC_{selected_run_type}_{beam_pass}_{selected_target_shortname}.csv")

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
    df["xsec_model"] = df["modelxsec"]
    df["bc_corr"] = 0.0

    df_out = df[["exp", "A", "Z","ebeam","theta","eprime","xbj", "q2", "w2", "epsilon","xsec_exp", "xsec_exp_err", "xsec_model", "bc_corr"]]

    df_out = df_out.dropna()

    all_rows.append(df_out)

df_data = pd.concat(all_rows, ignore_index=True)

output_csv = f"CSVs/{selected_run_type.upper()}_{selected_target_shortname}_xsec_results_compiled.csv"

# -----------------------------------------------------
# Starting logic to filter, figure out new binning
# -----------------------------------------------------
selected_var = "xbj"

df_xbj_minmax = (df_data.groupby("exp")[selected_var].agg(["min", "max"]).reset_index())
df_xbj_minmax = df_xbj_minmax.rename(columns={"min": f"{selected_var}_min","max": f"{selected_var}_max"})

overlap_min = df_xbj_minmax[f"{selected_var}_min"].max()
overlap_max = df_xbj_minmax[f"{selected_var}_max"].min()

overlap_range = float(overlap_max) - float(overlap_min)

stepsize = float(overlap_range) / float(nbins)

central_value = overlap_min + overlap_range / 2

if nbins % 2 == 1:
    #Odd: center is exact middle
    central_value = (overlap_min + overlap_max) / 2
    half_n = nbins // 2
    bin_centers = central_value + stepsize * np.arange(-half_n, half_n + 1)
else:
    #Even: two central bins around the middle
    central_value = (overlap_min + overlap_max) / 2
    half_n = nbins // 2
    bin_centers = central_value + stepsize * (np.arange(-half_n, half_n) + 0.5)

bin_centers = np.array(sorted(bin_centers))

edges = np.zeros(len(bin_centers)+1)

edges[1:-1] = (bin_centers[:-1] + bin_centers[1:]) / 2.0 #These are the center bins, they should be evenly spaced

edges[0]  = overlap_min  - (bin_centers[1]  - bin_centers[0]) / 2.0
edges[-1] = overlap_max + (bin_centers[-1] - bin_centers[-2]) / 2.0

centers_arr = bin_centers.reshape(1, -1)
vals = df_data[selected_var].to_numpy().reshape(-1, 1)
df_data["bin_num"] = np.argmin(np.abs(vals - centers_arr), axis=1) #I'm catching the 'lost' data points here: things outside the overlap range get assigned to closest bin.

print(f"Overlap range: {overlap_min:.4f} → {overlap_max:.4f}")
print(f"Stepsize: {stepsize:.4f}")
print(f"Bin centers ({nbins} bins):\n{bin_centers}")

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
        
    infile_name = f"{selected_run_type}_{beam_pass}_{selected_target_shortname}"
    
    for bin_center in bin_centers:
        input_string = f"1,{theta_inp:.3f},{bin_center:.6f},1,1\n{infile_name}\n"
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
                        
            if model_xsec_value is np.nan:
                print(f"Failed to parse xsec for {beam_pass}, bin {bin_center}")
                
            model_results.append({"A": selected_target_A,
                                  "Z": selected_target_Z,
                                  "ebeam": ebeam,
                                  "bc_theta": theta_inp,
                                  "bc_xbj": bin_center,
                                  "bc_q2": model_q2_value,
                                  "xsec_model": model_xsec_value})
        
        except Exception as e:
            print(f"Error running calc_dis_xsec: {e}")

df_centers = pd.DataFrame(model_results)
center_to_bin = dict(zip(bin_centers, range(len(bin_centers))))
df_centers["bin_num"] = df_centers["bc_xbj"].map(center_to_bin)

print(f"\nXSEC Model at Bin Centers")
print(df_centers)

df_centers_bin = (df_centers.groupby(["bin_num", "ebeam"], as_index=False)[["xsec_model", "bc_xbj", "bc_q2"]].mean())

df_data["bc_xsec_model"] = np.nan
df_data["bc_xbj"] = np.nan
df_data["bc_q2"] = np.nan

for _, row in df_centers_bin.iterrows():
    mask = (df_data["ebeam"] == row["ebeam"]) & (df_data["bin_num"] == row["bin_num"])
    df_data.loc[mask, "bc_xsec_model"] = row["xsec_model"]
    df_data.loc[mask, "bc_xbj"] = row["bc_xbj"]
    df_data.loc[mask, "bc_q2"] = row["bc_q2"]

df_data["bc_corr"] = df_data["bc_xsec_model"] / df_data["xsec_exp"]

col_final = ["exp", "A", "Z", "ebeam", "theta", "eprime", "xbj",
             "q2", "w2", "epsilon", "xsec_exp", "xsec_exp_err",
             "xsec_model", "bc_xsec_model", "bc_corr", "bin_num",
             "bc_xbj", "bc_q2"]

df_final = df_data[col_final]

df_final.to_csv(output_csv, index=False)

print(f"Saved compiled CSV → {output_csv}")
print(f"bc_corr range: {df_final['bc_corr'].min():.6f} → {df_final['bc_corr'].max():.6f}")

# -----------------------------------------------------
# Starting logic to filter, figure out new binning
# -----------------------------------------------------
plt.figure()

unique_bins = sorted(df_data["bin_num"].unique())

plt.axvspan(overlap_min, overlap_max, color="lightskyblue", alpha=0.1, label="Overlap Region")

plt.scatter(df_data["bc_xbj"],df_data["bc_q2"], marker="*",label="Bin centers", s=64)

for bin_num in unique_bins:
    mask = df_data["bin_num"] == bin_num
    plt.scatter(df_data.loc[mask, "xbj"],df_data.loc[mask, "q2"], label=f"Bin {bin_num}", s=12)

for i in range(1, len(edges)-1):
    if i == 1:
        plt.axvline(edges[i], ls="--", color="red", alpha=0.8, label="Bin Edges")
    else:
        plt.axvline(edges[i], ls="--", color="red", alpha=0.8)

# for i, value in enumerate([overlap_min, overlap_max]):
#     if i == 0:
#         plt.axvline(value, ls=":", color="red", alpha=0.4, label="Overlap Bounds")
#     else:
#         plt.axvline(value, ls=":", color="red", alpha=0.4)
            
    
# plt.axvline(overlap_min, ls=":", color="red", alpha=0.2, label="Overlap min/max")
# plt.axvline(overlap_max, ls=":", color="red", alpha=0.2)

plt.xlabel(r"x$_{bj}$")
plt.ylabel(r"Q$^2$")
plt.title(f"{selected_target_titlename} Binning Test (nbins={nbins})")

plt.legend()
plt.grid(axis="both", linestyle="--", alpha=0.8)

plt.savefig("PNGs/{selected_run_type)_{selected_target_shortname}_binning_test.png")

plt.show()

plt.close()
    



