#!/usr/bin/env python3

import os, re, csv, subprocess
from tqdm import tqdm
from datetime import datetime
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed

script_path = "HCAL_calib_check_skim.py"
bigtable_filepath = "../../../hallc_replay_rsidis/AUX_FILES/rsidis_bigtable_pass0.csv"
fit_results_filepath = "DAT/FIT_hcal_results.dat"
errlog_filepath = "ERR_log/FIT_hcal_error_log.txt"

def process_run(runnum, script_path):
    result = subprocess.run(["python3", script_path, runnum], capture_output=True, text=True, check=False, timeout=300)
    if result.returncode != 0:
        return None  # Silent fail, no print()
    fit_line = next((line.strip() for line in result.stdout.splitlines() if line.strip().startswith(runnum + "\t")), None)
    if fit_line is None:
        return None
    (run, mean, meanerr, sigma, sigmaerr, bin_min, bin_max, bin_avg, bin_sum) = fit_line.split("\t")

    return {run: [mean, meanerr, sigma, sigmaerr, bin_min, bin_max, bin_avg, bin_sum]}

if __name__ == '__main__':
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    
    runnums = []
    runinfo_lookup = {}
    if os.path.exists(bigtable_filepath):
        with open(bigtable_filepath, "r") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                try:
                    runnum = int(row["run"])
                    runtype = row.get("run_type", "").strip().lower()
                    runinfo_lookup[runnum] = {"runtype": runtype}
                    runnums.append(str(runnum))
                except (KeyError, ValueError):
                    continue
    else:
        raise FileNotFoundError(f"Missing bigtable: {bigtable_filepath}")

    updated_fit_results = {}
    errors = []
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(process_run, runnum, script_path): runnum for runnum in runnums}
        for future in tqdm(as_completed(futures), total=len(runnums)):
            runnum = futures[future]
            try:
                result = future.result()
                if result:
                    updated_fit_results.update(result)
            except Exception as e:
                errors.append((runnum, str(e)))

    # Error log
    with open(errlog_filepath, "w") as outfile:
        for runnum, error in errors:
            outfile.write(f"{runnum}: {error}\n")

    # Results
    with open(fit_results_filepath, "w") as outfile:
        outfile.write("runnum\truntype\tfit_mean\tmean_err\tfit_sigma\tsigma_err\t""bin_min\tbin_max\tbin_avg\tbin_total\n")

        for run in sorted(updated_fit_results, key=int):
            mean, meanerr, sigma, sigmaerr, bin_min, bin_max, bin_avg, bin_sum = updated_fit_results[run]
            runtype = runinfo_lookup.get(int(run), {}).get("runtype", "")
            outfile.write(f"{run}\t{runtype}\t{mean}\t{meanerr}\t{sigma}\t{sigmaerr}\t{bin_min}\t{bin_max}\t{bin_avg}\t{bin_sum}\n")

    tqdm.write(f"Success: {len(updated_fit_results)}/{len(runnums)}")
