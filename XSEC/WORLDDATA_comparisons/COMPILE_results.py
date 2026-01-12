#!/usr/bin/env python3

import pandas as pd
import os, re, sys

ratio_directory = "../FORM_xsec/RATIOS"

csv_files = [f"{ratio_directory}/XSEC_RATIO_hmsdis_4pass_c_to_ld2.csv",
             f"{ratio_directory}/XSEC_RATIO_hmsdis_5pass_c_to_ld2.csv",
             f"{ratio_directory}/XSEC_RATIO_hmsdis_4pass_cu_to_ld2.csv",
             f"{ratio_directory}/XSEC_RATIO_hmsdis_5pass_cu_to_ld2.csv"]

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

    df["Exp"] = f"RSIDIS {pass_label}Pass"
    df["A"] = df["num_A"]
    df["Z"] = df["num_Z"]
    df["x"] = df["xbj"]
    df["Q2"] = df["q2"]
    df["W2"] = df["w"] ** 2
    df["beamE(GeV)"] = df["eprime"]
    df["angle"] = df["theta"]
    df["Ratio"] = df["xsec_ratio_norm"]
    df["TotalErrorAbs"] = df["xsec_ratio_norm_err"]

    df_out = df[["Exp", "A", "Z", "x", "Q2", "W2", "beamE(GeV)", "angle", "Ratio", "TotalErrorAbs"]]

    df_out = df_out.dropna()

    all_rows.append(df_out)

df_all = pd.concat(all_rows, ignore_index=True)

output_csv = "RME_results.csv"

df_all.to_csv(output_csv, index=False)

print(f"Saved compiled CSV → {output_csv}")

    

