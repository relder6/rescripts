#!/usr/bin/env python3

import subprocess
from tqdm import tqdm
import os
import sys

script_path = "CALIB_CHECKS_pcal.py"
fit_results = []
fit_results_filepath = "FIT_results/FIT_results.dat"

# --------------------------------------------------------------------------
# Checking that fit results file exists
# --------------------------------------------------------------------------
if not os.path.exists(fit_results_filepath):
    print(f"ERROR: Output file {fit_results_filepath} does not exist! Exiting.")
    sys.exit(1)

# --------------------------------------------------------------------------
# Now looping over the runs
# --------------------------------------------------------------------------

# Full range of rsidis phase 1 runs is 23833, 25604
for runnum in tqdm(range(23833, 25604)):
    try:
        result = subprocess.run(
            ["python3", script_path],
            input=str(runnum),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text = True
            )
        lines = [line for line in result.stdout.splitlines() if line.strip()]

        fit_line = None
        for line in lines:
            if line.startswith(str(runnum)):
                fit_line = line
                break
        if fit_line is None:
            tqdm.write(f"WARNING: Fit results not found for run {runnum}. Skipping...")
            continue
        
        run, mean, sigma = fit_line.split("\t")
        fit_results.append((run, mean, sigma))
        # print(result.stdout.decode()) #Uncomment this line if you wish to see stdout prints.

    except subprocess.CalledProcessError as e:
        tqdm.write(f"ERROR code {e.returncode} for run {runnum}; check if this is an SHMS run.")
        continue

if not os.path.exists(fit_results_filepath):
    print(f"ERROR: Output file {fit_results_filepath} does not exist! Exiting.")
    sys.exit(1)

# --------------------------------------------------------------------------
# Opening existing stats file
# --------------------------------------------------------------------------
    
with open(fit_results_filepath, "r") as infile:
    lines = infile.readlines()

header = lines[0]
fit_result_lines = {line.split("\t")[0]: line.strip().split("\t")[1:] for line in lines[1:] if line.strip()}

for run, mean, sigma in fit_results:
    fit_result_lines[run] = [mean, sigma]

# --------------------------------------------------------------------------
# Writes fit results to stats file
# --------------------------------------------------------------------------
with open(fit_results_filepath, "w") as outfile:
    outfile.write(header)
    for run, values in sorted(fit_result_lines.items(), key=lambda x: int(x[0])):
        outfile.write(f"{run}\t{values[0]}\t{values[1]}\n")
