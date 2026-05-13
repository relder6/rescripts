#!/usr/bin/env python3

import os, re, sys, csv
import numpy as np
import subprocess
import pandas as pd
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import get_data_cuts, get_common_values
from INIT.config import parse_run_type, parse_beam_pass, parse_target, parse_bins

# -----------------------------------------------------
# Handling user inputs, listing directories
# -----------------------------------------------------
arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None
arg3 = sys.argv[3] if len(sys.argv) > 3 else None

selected_run_type = parse_run_type(arg1)
selected_beam_pass, beam_prefix = parse_beam_pass(arg2)
target_abbrev, target_longname, target_shortname, target_A, target_Z = parse_target(arg3)

vals = get_common_values()
ebeam_4pass = vals["ebeam_4pass"]
theta_4pass = vals["angle_4pass"]
ebeam_5pass = vals["ebeam_5pass"]
theta_5pass = vals["angle_5pass"]

model_xsec_dir = "../../../mc-single-arm/util/dis_xec"

infile_name = f"{selected_run_type}_{selected_beam_pass}pass_{target_abbrev}"

outfile_name = f"{selected_run_type}_{selected_beam_pass}pass_{target_abbrev}_model_xsec.csv"
# outfile_name = "testing.csv"
# -----------------------------------------------------
# Now building input strings, collecting model xsec of bin centers
# -----------------------------------------------------
model_results = []
found_header = False

if selected_beam_pass == "4":
    theta_inp = theta_4pass
    ebeam = ebeam_4pass
    ep_min = 1.385555
    ep_step = 0.01531
    ep_num = 20
    p0 = 1.531
    
elif selected_beam_pass == "5":
    theta_inp = theta_5pass
    ebeam = ebeam_5pass
    ep_min = 3.296010
    ep_step = 0.036420
    ep_num = 20
    p0 = 3.642

theta_rad = np.radians(float(theta_inp))

# The input string depends on the version of Dave's xsec tool, mc-single-arm/util/dis_xec/calc_dis_xsec
# Right now, the input string is flag (0 = fixed theta, bin in eprime; 1 = fixed theta, bin in xbj; 2 = fixed Q2, bin in xbj...
# fixed var (theta or Q2), <var>min, <var>step, <var>bin_num, then input filename
input_string = f"0\n{theta_inp:.6f}\n{ep_min:.6f}\n{ep_step:.6f}\n{ep_num}\n{infile_name}\n"
try:
    result = subprocess.run(["./calc_dis_xsec"], input=input_string, text=True, capture_output=True, cwd=model_xsec_dir)
    print("STDOUT:")
    print(result.stdout)
    print("STDERR:")
    print(result.stderr)
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
        if not found_header:
            continue

        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            model_eprime = float(parts[0])
            model_theta = float(parts[1])
            model_xbj = float(parts[2])
            model_q2 = float(parts[3])
            model_w = float(parts[4])
            model_xsec = float(parts[-1])
        except ValueError:
            continue
        nu = ebeam - model_eprime

        if model_q2 > 0:
            epsilon = 1 / (1 + 2 * (nu**2 + model_q2) / model_q2 * np.tan(theta_rad/2)**2)
        else:
            epsilon = np.nan

        delta = round(100 * (model_eprime - p0) / p0, 2)
        # eprime,theta,xbj,q2,w,modelxsec,epsilon,delta    
        model_results.append({"eprime": model_eprime,
                              "theta": model_theta,
                              "xbj": model_xbj,
                              "q2": model_q2,
                              "w": model_w,
                              "modelxsec": model_xsec,
                              "epsilon": epsilon,
                              "delta": delta})
        
except Exception as e:
    print(f"Error running calc_dis_xsec: {e}")

df = pd.DataFrame(model_results)

with open(outfile_name, "w", newline="") as outfile:
    writer = csv.DictWriter(outfile, fieldnames = model_results[0].keys())

    writer.writeheader()
    writer.writerows(model_results)

print(f"Wrote {outfile_name}")
