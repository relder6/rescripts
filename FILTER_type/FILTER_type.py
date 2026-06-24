#!/usr/bin/env python3

import os, re, sys, csv
import numpy as np

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import parse_run_type, parse_beam_pass, parse_target, get_flags

# -----------------------------------------------------
# Handling user inputs
# -----------------------------------------------------
flags = get_flags()

USING_CURRENT_OFFSET = flags["USING_CURRENT_OFFSET"]
USING_BOIL_CORR = flags["USING_BOIL_CORR"]
USING_CURRENT_CUT = flags["USING_CURRENT_CUT"]

if USING_CURRENT_CUT:
    current_cut = 10.0
else:
    current_cut = 0.0

skip_runnums = [23853, 23854, 23855, 23856, 23857, 23858, 23859, 23860]
                #Now skipping the low current (< 10 uA) runs,
                # 23934,23938,23963,24027,24290,24291,24292,24293,24294,24308,24309,
                # 24319,24333,24438,24440,24455,24456,24481,24482,24483,24495,24496,
                # 24498,24911,24967,25047,25081,25406,25407,25416,25417,
                #Now skipping some runs that have negative yields (?!),
                # 25396, 25397]

arg1 = sys.argv[1] if len(sys.argv) > 1 else None
arg2 = sys.argv[2] if len(sys.argv) > 2 else None
arg3 = sys.argv[3] if len(sys.argv) > 3 else None

selected_run_type = parse_run_type(arg1)
selected_beam_pass, beam_prefix = parse_beam_pass(arg2)
target_abbrev, target_longname, target_shortname, target_A, target_Z = parse_target(arg3)

bigtable_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_bigtable_pass0p1.csv"
output_filepath = f"{target_abbrev.upper()}/{selected_run_type}_{selected_beam_pass}pass_{target_abbrev}_runs.csv"

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
                ibeam1 = row.get("BCM1_I", "N/A")
                ibeam4a = row.get("BCM4A_I", "N/A")
                ibeam4c = row.get("BCM4C_I", "N/A")
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
                if target_abbrev in ["lh2", "ld2"]:
                    if USING_BOIL_CORR:
                        boil_corr = row.get("boil_corr", "N/A")
                    else:
                        boil_corr = 1.0
                else:
                    boil_corr = 1.0

                # Calculating extra values (like weights) that downstream scripts expect,
                if (selected_run_type == run_type.strip().lower() and
                    target_abbrev.lower() == target.strip().lower() and
                    ebeam.startswith(beam_prefix) and
                    runnum not in skip_runnums):

                    if ibeam in ["N/A", ""]:
                        continue
                    try:
                        if float(ibeam) < current_cut:
                            continue
                    except (ValueError, TypeError):
                        continue

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

                    if USING_CURRENT_OFFSET:
                        current_offset = -0.0301
                        current_offset_corr = 1 / (1 + (current_offset / float(ibeam)))
                    else:
                        current_offset_corr = 1

                    weight = (float(boil_corr) * float(ps) * float(current_offset_corr)) / (float(livetime) * float(trackeff))

                    if runnum <= 27100:
                        phase = "I"
                    elif runnum > 27100:
                        phase = "II"

                    bigtable_lookup[runnum] = {"run_type": run_type,
                                               "start_time": start_time,
                                               "stop_time": stop_time,
                                               "boil_corr": boil_corr,
                                               "fan_mean": fan_mean,
                                               "ebeam": ebeam,
                                               "ibeam": ibeam,
                                               "ibeam1": ibeam1,
                                               "ibeam4a": ibeam4a,
                                               "ibeam4c": ibeam4c,
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
                                               "weight": weight,
                                               "phase": phase}
            except (ValueError, KeyError):
                continue
                               
# -----------------------------------------------------
# Writing output table
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
                target_abbrev.lower() == data["target"].strip().lower() and
                data["ebeam"].startswith(beam_prefix) and
                runnum not in skip_runnums):
            
                writer.writerow({"runnum": runnum, **data})
                rows_written += 1

    print(f"\n☢️  Saved {rows_written} matching lines to {output_filepath}")

else:
    print(f"\n⚠️  No runs matched '{selected_run_type},{selected_beam_pass}Pass,{target_abbrev}'. No file was written.")
