#!/usr/bin/env python3

import os, re
import numpy as np
import csv

# =====================================================================
# Defining paths here
# =====================================================================

input_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat" # The location of the auxfiles runlist
hms_report_filepath = "/work/hallc/c-rsidis/replay/pass0/REPORT_OUTPUT/HMS/PRODUCTION/replay_hms_coin_production_{runnum}_-1.report"
shms_report_filepath = "/work/hallc/c-rsidis/replay/pass0/REPORT_OUTPUT/SHMS/PRODUCTION/replay_shms_coin_production_{runnum}_-1.report"
coin_report_filepath = "/work/hallc/c-rsidis/replay/pass0/REPORT_OUTPUT/COIN/PRODUCTION/replay_coin_production_{runnum}_-1.report"
run_info_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_bigtable_pass0.csv"

skip_runnums = [23853, 23854, 23855, 23856, 23857, 23858, 23859, 23860, 24482, 24496, 25013]
# Runs 23853 - 23860 are commissioning runs; skipping these for now
# What is happening with yields on 24482, 24496, 25013?

# =====================================================================
# Handling user inputs
# =====================================================================

selected_type = input("Enter desired run type (default HMSDIS): ").strip().lower()
if not selected_type:
    selected_type = f"hmsdis"
    
selected_beam_pass = input("Enter desired beam pass (present options: 1, 4, 5): ").strip()
selected_beam_pass_to_energy_prefix = {"1": "2.",
                                       "2": "4.",
                                       "3": "6.",
                                       "4": "8.",
                                       "5": "10."}
beam_prefix = selected_beam_pass_to_energy_prefix.get(selected_beam_pass)
if not beam_prefix:
    print(f"Unknown pass: {selected_beam_pass}.  Please try again.")
    exit(1)

# =====================================================================
# Reading auxfiles runlist, filtering and extracting lines
# =====================================================================

with open(input_filepath, "r") as infile:
    lines = infile.readlines()

run_lines = []

for line in lines:
    if line.lstrip().startswith("#"):        # dropping comment lines
        continue
    if line.lstrip().startswith("!"):        # dropping lines starting with !
        continue
    if re.match(r"^\s*[-=*]", line):         # dropping separator lines
        continue
    run_lines.append(line)

filtered_lines = []

for line in run_lines:
    parts = re.split(r'\s+', line.strip())
    ebeam, target_type, run_type = "", "", ""
    for part in parts:
        if re.match(r'^\d+\.\d+$', part):
            ebeam = part
            break
    for part in parts:
        if part.lower() in ["hole", "optics", "heep", "hee", "hmsdis", "shmsdis", "pi-sidis", "pi+sidis", "junk"]:
            run_type = part.lower()
            break
    for part in parts:
        if part.lower() in ["lh2", "c", "cu", "al", "ld2", "dummy", "hole"]:
            target_type = part.lower()
            break
    runnum = int(parts[0])
    beam_match = ebeam.startswith(beam_prefix)

    if  beam_match and runnum not in skip_runnums:
        filtered_lines.append(line)

output_filepath = f"{selected_beam_pass}pass_rates.tsv"

# =====================================================================
# Reading in report files and compiling the tsv
# =====================================================================
        
with open(output_filepath, "w") as outfile:
    outfile.write("#Run#\tRunType\tibeam\tHMS3of4\tSHMS3of4\n")
    
    for line in filtered_lines:
        parts = re.split(r'\s+', line.strip())
        runtype = parts[11].lower()
        runnum = parts[0]
        ibeam = parts[4]

        hms3of4 = -999
        shms3of4 = -999
        
        if runtype == 'hmsdis':
            hms_report_file = hms_report_filepath.format(runnum=runnum)
            if os.path.exists(hms_report_file):
                with open(hms_report_file, "r") as rep:
                    for rep_line in rep:
                        if rep_line.startswith("HMS 3/4 Trigger Rate"):
                            m = re.search(r":\s*([0-9.]+)", rep_line)
                            if m:
                                hms3of4 = float(m.group(1))

        if runtype == 'shmsdis':
            shms_report_file = shms_report_filepath.format(runnum=runnum)
            if os.path.exists(shms_report_file):
                with open(shms_report_file, "r") as rep:
                    for rep_line in rep:
                        if rep_line.startswith("SHMS 3/4 Trigger Rate"):
                            m = re.search(r":\s*([0-9.]+)", rep_line)
                            if m:
                                shms3of4 = float(m.group(1))
                                
        if runtype in ('pi+sidis','pi-sidis'):
            coin_report_file = coin_report_filepath.format(runnum=runnum)
            if os.path.exists(coin_report_file):
                with open(coin_report_file, "r") as rep:
                    for rep_line in rep:
                        
                        if rep_line.startswith("SHMS 3/4 Trigger Rate"):
                            m = re.search(r":\s*([0-9.]+)", rep_line)
                            if m:
                                shms3of4 = float(m.group(1))

                        if rep_line.startswith("HMS 3/4 Trigger Rate"):
                            m = re.search(r":\s*([0-9.]+)", rep_line)
                            if m:
                                hms3of4 = float(m.group(1))
        tsv_line = "\t".join(map(str, [runnum, runtype, ibeam, f"{hms3of4:.8g}", f"{shms3of4:.8g}"]))
        outfile.write(tsv_line + "\n")
            

print(f"\n☢️  Wrote {len(filtered_lines)} matching lines to {output_filepath}")
