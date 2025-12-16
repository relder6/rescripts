#!/usr/bin/env python3

import os, re, sys
import numpy as np

input_file = sys.argv[1] if len(sys.argv) > 1 else input("Input parameter file: ")
input_filename = f"{input_file}.dat"
output_filename = f"{input_file}.csv"

if not os.path.exists(input_filename):
    print(f"ERROR: '{input_filename}' not found.  Exiting...")
    sys.exit(1)

with open(input_filename, "r") as infile:
    with open(output_filename, "w") as outfile:
        outfile.write("eprime,theta,xbj,q2,w,modelxsec\n")
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

            outfile.write(",".join([f"{c}" if not np.isnan(c) else "nan" for c in cols]) + "\n")
            
print(f"Made {output_filename}.")

