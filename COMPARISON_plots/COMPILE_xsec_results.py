#!/usr/bin/env python3

import pandas as pd
import os, re, sys

ratio_directory = "../XSEC/FORM_xsec"

csv_files = [f"{ratio_directory}/C/XSEC_hmsdis_4pass_c.csv",
             f"{ratio_directory}/C/XSEC_hmsdis_5pass_c.csv",
             f"{ratio_directory}/CU/XSEC_hmsdis_4pass_cu.csv",
             f"{ratio_directory}/CU/XSEC_hmsdis_5pass_cu.csv",
             f"{ratio_directory}/LD2/XSEC_hmsdis_4pass_ld2.csv",
             f"{ratio_directory}/LD2/XSEC_hmsdis_5pass_ld2.csv",
             f"{ratio_directory}/LH2/XSEC_hmsdis_4pass_lh2.csv",
             f"{ratio_directory}/LH2/XSEC_hmsdis_5pass_lh2.csv"]

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
        ebeam = 8.5831
    elif int(pass_label) == 5:
        ebeam = 10.6716

    # A,Z,eprime,theta,xbj,q2,w,epsilon,modelxsec,xsec_exp,xsec_exp_err

    df["exp"] = f"RSIDIS({pass_label}Pass)"
    df["ebeam"] = ebeam
    df["q2"] = df["q2"]
    df["w2"] = df["w"] ** 2
    df["xsec_model"] = df["modelxsec"]

    df_out = df[["exp", "A", "Z","ebeam","eprime","theta","xbj", "q2", "w2", "epsilon","xsec_exp", "xsec_exp_err", "xsec_model"]]

    df_out = df_out.dropna()

    all_rows.append(df_out)

df_all = pd.concat(all_rows, ignore_index=True)

output_csv = "CSVs/RME_xsec_results.csv"

df_all.to_csv(output_csv, index=False)

print(f"Saved compiled CSV → {output_csv}")

    

