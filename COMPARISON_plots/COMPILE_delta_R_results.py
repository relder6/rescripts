#!/usr/bin/env python3

import pandas as pd
import os, re, sys

R_directory = "../LT_separations/DELTA_R/CSVs"

csv_files = [f"{R_directory}/DELTA_R_HMSDIS_rosenbluth_fit_c_to_ld2.csv",
             f"{R_directory}/DELTA_R_HMSDIS_rosenbluth_fit_cu_to_ld2.csv",
             f"{R_directory}/DELTA_R_HMSDIS_rosenbluth_fit_al_to_ld2.csv"]

all_rows = []

# bin_num,xbj,q2_avg,q2_err,delta_R,delta_R_err

for filepath in csv_files:
    if not os.path.exists(filepath):
        print(f"Skipping missing file: {filepath}")
        continue

    df = pd.read_csv(filepath)

    m = re.search(r'fit_([a-z0-9]+)\_to_ld2.csv', os.path.basename(filepath))
    if m:
        target_label = m.group(1)
    else:
        print(f"Could not determine target from filename: {filepath}, defaulting to unknown.")
        target_label = "???"

    df["exp"] = f"RSIDIS"
    df["target"] = target_label
    df["bin_num"] = df["bin_num"]
    
    df_out = df[["exp", "target", "bin_num", "xbj", "q2_avg", "delta_R", "delta_R_err", "sigma_t_ratio", "sigma_t_ratio_err"]]

    df_out = df_out.dropna()

    all_rows.append(df_out)

df_all = pd.concat(all_rows, ignore_index=True)

output_csv = "CSVs/RME_Delta_R_results.csv"

df_all.to_csv(output_csv, index=False)

print(f"Saved compiled CSV → {output_csv}")

    

