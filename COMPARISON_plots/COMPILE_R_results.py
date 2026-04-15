#!/usr/bin/env python3

import pandas as pd
import os, re, sys

R_directory = "../LT_separations/NUCLEAR_R/CSVs"

csv_files = [f"{R_directory}/HMSDIS_rosenbluth_fit_c.csv",
             f"{R_directory}/HMSDIS_rosenbluth_fit_cu.csv",
             f"{R_directory}/HMSDIS_rosenbluth_fit_ld2.csv",
             f"{R_directory}/HMSDIS_rosenbluth_fit_lh2.csv"]

all_rows = []

# bin_num,xbj,q2_avg,q2_err,sigma_T,sigma_T_err,sigma_L,sigma_L_err,R,R_err,n_eps

for filepath in csv_files:
    if not os.path.exists(filepath):
        print(f"Skipping missing file: {filepath}")
        continue

    df = pd.read_csv(filepath)

    m = re.search(r'fit_([a-z0-9]+)\.csv', os.path.basename(filepath))
    if m:
        target_label = m.group(1)
    else:
        print(f"Could not determine target from filename: {filepath}, defaulting to unknown.")
        target_label = "???"

    df["exp"] = f"RSIDIS"
    df["target"] = target_label
    df["bin_num"] = df["bin_num"]
    
    df_out = df[["exp", "target", "bin_num", "R", "R_err", "n_eps"]]

    df_out = df_out.dropna()

    all_rows.append(df_out)

df_all = pd.concat(all_rows, ignore_index=True)

output_csv = "CSVs/RME_Nuclear_R_results.csv"

df_all.to_csv(output_csv, index=False)

print(f"Saved compiled CSV → {output_csv}")

    

