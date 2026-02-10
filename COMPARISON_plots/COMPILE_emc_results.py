#!/usr/bin/env python3

import pandas as pd
import os, re, sys

ratio_directory = "../XSEC/FORM_xsec/RATIOS"

csv_files = [f"{ratio_directory}/XSEC_RATIO_hmsdis_4pass_c_to_ld2.csv",
             f"{ratio_directory}/XSEC_RATIO_hmsdis_5pass_c_to_ld2.csv",
             f"{ratio_directory}/XSEC_RATIO_hmsdis_4pass_cu_to_ld2.csv",
             f"{ratio_directory}/XSEC_RATIO_hmsdis_5pass_cu_to_ld2.csv"]

all_rows = []

# A_num,Z_num,eprime,theta,xbj,q2,w,epsilon,modelxsec_num,xsec_exp_num,xsec_exp_err_num,A_denom,Z_denom,modelxsec_denom,xsec_exp_denom,xsec_exp_err_denom,xsec_ratio,xsec_ratio_err,xsec_ratio_norm,xsec_ratio_norm_err,f2rat,iso_corr,xsec_ratio_final,xsec_ratio_final_err

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

    df["exp"] = f"RSIDIS({pass_label}Pass)"
    df["A"] = df["A_num"]
    df["Z"] = df["Z_num"]
    df["w2"] = df["w"] ** 2
    df["ratio_raw"] = df["xsec_ratio"]
    df["ratio_norm"] = df["xsec_ratio_norm"]
    df["emc_ratio"] = df["xsec_ratio_final"]
    df["emc_ratio_err"] = df["xsec_ratio_final_err"]
    

    df_out = df[["exp", "A", "Z", "xbj", "q2", "w2", "epsilon", "eprime", "theta", "xsec_exp_num", "xsec_exp_denom", "ratio_raw", "ratio_norm", "iso_corr", "emc_ratio", "emc_ratio_err"]]

    df_out = df_out.dropna()

    all_rows.append(df_out)

df_all = pd.concat(all_rows, ignore_index=True)

output_csv = "CSVs/RME_emc_results.csv"

df_all.to_csv(output_csv, index=False)

print(f"Saved compiled CSV → {output_csv}")

    

