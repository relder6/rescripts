#!/usr/bin/env python3

import os, re, sys, csv
import numpy as np

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import get_common_run_inputs

# -----------------------------------------------------
# Establishing paths
# -----------------------------------------------------
input_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat" # The location of the auxfiles runlist
report_filepath = "/work/hallc/c-rsidis/replay/pass0/REPORT_OUTPUT/HMS/PRODUCTION/replay_hms_coin_production_{runnum}_-1.report"
bigtable_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_bigtable_pass0.csv"

skip_runnums = [23853, 23854, 23855, 23856, 23857, 23858, 23859, 23860,
                25396, 25397,
                23918,23919,23934,23938,23963,24027,24290,24291,24292,24293,24294,24308,24309,
                24319,24333,24438,24440,24455,24456,24481,24482,24483,24495,24496,24498, 
                24911,24967,25047,25081,25406,25407,25416,25417]
# Runs 23853 - 23860 are commissioning runs; skipping these for now
# Runs 25396 and 25397 have NEGATIVE YIELDS in my scripts for some reason???!!!
# Low current runs, 4pass:
# 23918,23919,23934,23938,23963,24027,24290,24291,24292,24293,24294,24308,24309,24319,24333,24438,24440,24455,24456,24481,24482,24483,24495,24496,24498
# Low current runs, 5pass:
# 24911,24967,25047,25081,25406,25407,25416,25417

# -----------------------------------------------------
# Input handling
# -----------------------------------------------------
selected_run_type, selected_beam_pass, beam_prefix, selected_target_shortname, selected_target_titlename, selected_target_A, selected_target_Z = get_common_run_inputs()

output_filepath = f"{selected_target_shortname.upper()}/{selected_run_type}_{selected_beam_pass}pass_{selected_target_shortname}_runs.dat"
# output_filepath = "testing.csv"
# -----------------------------------------------------
# Bigtable look-up
# -----------------------------------------------------
bigtable_lookup = {}
if os.path.exists(bigtable_filepath):
    with open(bigtable_filepath, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                # Purely reading from the bigtable here,
                runnum = int(row["run"])
                run_type = row.get("run_type", "N/A")
                target = row.get("target", "N/A")
                ebeam = row.get("ebeam", "N/A")
                ibeam = row.get("BCM2_I", "N/A")
                qbeam = row.get("BCM2_Q", "N/A")
                hms_p = row.get("hms_p", "N/A")
                hms_th = row.get("hms_th", "N/A")
                ps3 = row.get("ps3", "N/A")
                ps4 = row.get("ps4", "N/A")
                trackeff = row.get("h_esing_Eff", "N/A")
                livetime = row.get("comp_livetime", "N/A")
                fan_mean = row.get("fan_mean", "N/A")
                start_time = row.get("start_time", "N/A")
                stop_time = row.get("stop_time", "N/A")

                # Being sure to only read out the boil_corr for liquid targets,
                if selected_target_shortname in ["lh2", "ld2"]:
                    boil_corr = row.get("boil_corr", "N/A")
                else:
                    boil_corr = 1.0

                # Calculating extra values (like weights) that downstream scripts expect,
                if (selected_run_type == run_type.strip().lower() and
                    selected_target_shortname.lower() == target.strip().lower() and
                    ebeam.startswith(beam_prefix) and
                    runnum not in skip_runnums):

                    ps3_val, ps4_val = float(ps3), float(ps4)

                    if ps3_val == -999 or ps4_val == -999:
                        print(f"WARNING: run {runnum} has issues with prescale values.  Fix the lookup table for this setting!")
                        continue
                    elif ps3_val != -1 and ps4_val != -1:
                        print(f"WARNING: run {runnum} has both prescales set (ps3={ps3_val}, ps4={ps4_val}), using ps3.")
                    elif ps3_val != -1 and ps4_val == -1:
                        ps = ps3_val
                    elif ps4_val != -1 and ps3_val == -1:
                        ps = ps4_val
                    else:
                        print(f"Run {runnum} has no valid prescale (both -1), skipping...")
                        continue
                    
                   
                    nu = abs(float(ebeam)) - abs(float(hms_p))
                    hms_th_rad = np.deg2rad(float(hms_th))
                    q2 = 4 * abs(float((ebeam)) * abs(float(hms_p)) * (np.sin(hms_th_rad/2))**2)
                    epsilon = 1 / (1 + 2 * (1 + (nu**2 / q2)) * np.tan(hms_th_rad/ 2))**2
                    weight = float(boil_corr) * float(ps) / (float(livetime) * float(trackeff))

                    bigtable_lookup[runnum] = {"run_type": run_type,
                                               "start_time": start_time,
                                               "stop_time": stop_time,
                                               "boil_corr": boil_corr,
                                               "fan_mean": fan_mean,
                                               "ebeam": ebeam,
                                               "ibeam": ibeam,
                                               "qbeam": qbeam,
                                               "target": target,
                                               "hms_p": hms_p,
                                               "hms_th": hms_th,
                                               "hms_th_rad": hms_th_rad,
                                               "ps3": ps3,
                                               "ps4": ps4,
                                               "trackeff": trackeff,
                                               "nu": nu,
                                               "q2": q2,
                                               "epsilon": epsilon,
                                               "weight": weight}
            except (ValueError, KeyError):
                continue
                    
            
# -----------------------------------------------------
# Bigtable look-up
# -----------------------------------------------------
output_dir = os.path.dirname(output_filepath)
if output_dir:
    os.makedirs(output_dir, exist_ok=True)

if bigtable_lookup:
    columns = ["runnum"] + list(next(iter(bigtable_lookup.values())).keys())

    with open(output_filepath, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns)
        writer.writeheader()

        rows_written = 0

        for runnum, data in bigtable_lookup.items():

            if (selected_run_type == data["run_type"].strip().lower() and
                selected_target_shortname.lower() == data["target"].strip().lower() and
                data["ebeam"].startswith(beam_prefix) and
                runnum not in skip_runnums):
            
                writer.writerow({"runnum": runnum, **data})
                rows_written += 1

    print(f"\n☢️  Saved {rows_written} matching lines to {output_filepath}")

else:
    print(f"\n⚠️  No runs matched '{selected_run_type},{selected_beam_pass}Pass,{selected_target_shortname}'. No file was written.")
