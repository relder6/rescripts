#!/usr/bin/env python3

import os, re, csv, subprocess
from tqdm import tqdm

script_path = "CALIB_CHECKS_hcal.py"
fit_results_filepath = "DAT/FIT_hcal_results.dat"
auxfiles_runlist_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat"

# --------------------------------------------------------------------------
# Load existing fit results
# --------------------------------------------------------------------------
fit_result_lines = {}
if os.path.exists(fit_results_filepath):
    with open(fit_results_filepath, "r") as infile:
        lines = infile.readlines()
        if len(lines) > 1:
            for line in lines[1:]:
                if line.strip():
                    runnum, mean, meanerr, sigma, sigmaerr = line.strip().split("\t")
                    fit_result_lines[runnum] = [mean, meanerr, sigma, sigmaerr]
                    
# --------------------------------------------------------------------------
# Extracting valid run numbers from the auxfiles runlist
# --------------------------------------------------------------------------

def extract_parts(line):
    parts = re.split(r'\s+', line.strip())
    if not parts:
        return None
    if len(parts) >=12:
        run_type = parts[11].strip().lower()
        if run_type == "junk":
            return None
    return parts[0]

runnums = []

with open(auxfiles_runlist_filepath, "r") as infile:
    for line in infile:
        line = line.strip()
        if not line or line.startswith(("#", "!", "-", "=", "*")):
            continue
        runnum = extract_parts(line)
        if runnum:
            runnums.append(runnum)

# --------------------------------------------------------------------------
# Processing runs
# --------------------------------------------------------------------------
updated_fit_results = {}

for runnum in tqdm(runnums):
    try:
        result = subprocess.run(
            ["python3", script_path, runnum],
            capture_output=True, text=True, check=False)

        if result.returncode !=0:
            tqdm.write(f"WARNING: Run {runnum} failed (returncode {result.returncode}): skipping...")
            continue

        fit_line = next((line for line in result.stdout.splitlines() if line.startswith(runnum + "\t")), None)

        if fit_line is None:
            tqdm.write(f"WARNING: Fit results not found for run {runnum}; skipping...")
            continue

        run, mean, meanerr, sigma, sigmaerr = fit_line.split("\t")
        updated_fit_results[run] = [mean, meanerr, sigma, sigmaerr]

    except Exception as e:
        tqdm.write(f"ERROR: Unexpected error for run {runnum}: {e}")
        continue

# --------------------------------------------------------------------------
# Write fit results back to file
# --------------------------------------------------------------------------
with open(fit_results_filepath, "w") as outfile:
    outfile.write("#runnum\tfit_mean\tmean_err\tfit_sigma\tsigma_err\n")
    for run, (mean, meanerr, sigma, sigmaerr) in sorted(updated_fit_results.items(), key=lambda x: int(float(x[0]))):
        outfile.write(f"{run}\t{mean}\t{meanerr}\t{sigma}\t{sigmaerr}\n")

tqdm.write(f"Processing complete. Total runs processed: {len(updated_fit_results)}")
