#!/usr/bin/env python3

import os, re, sys
import numpy as np

input_file = sys.argv[1] if len(sys.argv) > 1 else input("Input parameter file: ")
input_filename = f"{input_file}.dat"
output_filename = f"{input_file}.csv"

if not os.path.exists(input_filename):
    print(f"ERROR: '{input_filename}' not found.  Exiting...")
    sys.exit(1)

# Determine beam energy based on filename
if "4pass" in input_file:
    E_beam = 8.5831  # GeV
    theta_deg = 29.045 #deg
elif "5pass" in input_file:
    E_beam = 10.6716  # example, replace with correct value if different
    theta_deg = 16.75 #deg
else:
    print(f"Beam pass not detected from input file name, exiting...")
    exit(1)
    
theta_rad = np.radians(theta_deg)

with open(input_filename, "r") as infile:
    with open(output_filename, "w") as outfile:
        outfile.write("eprime,theta,xbj,q2,w,modelxsec,epsilon\n")
        for line in infile:
            if line.lstrip().startswith("#") or line.lstrip().startswith(";"):
                continue
            cols = line.strip().split()
            if len(cols) != 6:
                cols = [np.nan]*6
            for i in range(6):
                try:
                    cols[i] = float(cols[i])
                except ValueError:
                    cols[i] = np.nan

            eprime = cols[0]
            q2 = cols[3]

            # Compute energy transfer
            nu = E_beam - eprime

            # Compute epsilon
            if not np.isnan(eprime) and not np.isnan(q2):
                epsilon = 1 / (1 + 2 * (nu**2 + q2)/q2 * np.tan(theta_rad/2)**2)
            else:
                epsilon = np.nan

            outfile.write(",".join([f"{c}" if not np.isnan(c) else "nan" for c in cols] + [f"{epsilon}"]) + "\n")
            
print(f"Wrote {output_filename}")
