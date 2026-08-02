#!/usr/bin/env python3

import os, re, sys, csv
import numpy as np
import pandas as pd
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)    
from INIT.config import parse_run_type, parse_beam_pass, parse_target, get_flags, parse_phase, parse_setting, get_common_values

flags = get_flags()
vals = get_common_values()

USING_CURRENT_OFFSET = flags["USING_CURRENT_OFFSET"]
USING_BOIL_CORR = flags["USING_BOIL_CORR"]
USING_CURRENT_CUT = flags["USING_CURRENT_CUT"]

# -----------------------------------------------------
# Things to skip, loop over, etc.
# -----------------------------------------------------
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
                
bigtable_filepaths = {"I": "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_bigtable_pass0p1.csv",
                      "II": "/lustre24/expphy/volatile/hallc/c-rsidis/relder/STUFF/rsidis_bigtable_phaseII.csv"}
report_filepaths = {"I": "/w/hallc-scshelf2102/c-rsidis/replay/pass0p1/REPORT_OUTPUT/HMS/PRODUCTION",
                    "II": "/volatile/hallc/c-rsidis/relder/STUFF/REPORT_PHASEII"}

master_output_filepath = f"MASTER_hmsdis.csv"

# -----------------------------------------------------
# Bigtable look-up
# -----------------------------------------------------
targets = ["al", "c", "cu", "dummy", "ld2", "lh2"]
angle_tolerance = 0.2
phases = ["I", "II"]
bigtable_rows_master = []

for phase in phases:
    bigtable_filepath = bigtable_filepaths[phase]
    report_filepath = report_filepaths[phase]

    if os.path.exists(bigtable_filepath):
        with open(bigtable_filepath, "r") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                try:
                    runnum = int(row["run"])
                    run_type = row.get("run_type", "N/A")
                    if run_type.strip().lower() != "hmsdis":
                        continue
                    target = row.get("target", "N/A")
                    ebeam = float(row.get("ebeam", "N/A"))
                    ibeam_1 = float(row.get("BCM1_I", "N/A"))
                    ibeam_2 = float(row.get("BCM2_I", "N/A"))
                    ibeam_4a = float(row.get("BCM4A_I", "N/A"))
                    ibeam_4c = float(row.get("BCM4C_I", "N/A"))
                    qbeam_1 = float(row.get("BCM1_Q", "N/A"))
                    qbeam_2 = float(row.get("BCM2_Q", "N/A"))
                    qbeam_4a = float(row.get("BCM4A_Q", "N/A"))
                    qbeam_4c = float(row.get("BCM4C_Q", "N/A"))
                    hms_p = float(row.get("hms_p", "N/A"))
                    hms_th = float(row.get("hms_th", "N/A"))
                    ps3 = row.get("ps3", "N/A")
                    ps4 = row.get("ps4", "N/A")
                    trackeff = float(row.get("h_esing_Eff", "N/A"))
                    livetime = float(row.get("comp_livetime", "N/A"))

                    # Being sure to only read out the boil_corr for liquid targets,
                    if target.lower() in ["lh2", "ld2"]:
                        if USING_BOIL_CORR:
                            boil_corr = row.get("boil_corr", "N/A")
                        else:
                            boil_corr = 1.0
                    else:
                        boil_corr = 1.0

                    try:
                        if float(ibeam_2) < current_cut:
                            continue
                    except (ValueError, TypeError):
                        continue
                    try:
                        ps3_val, ps4_val = float(ps3), float(ps4)
                    except (ValueError, TypeError):
                        print(f"WARNING: run {runnum} has non-numeric prescales (ps3 = {ps3}, ps4 = {ps4}), skipping...")
                        continue

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

                    mp = 0.93827208943 #in GeV
                    nu = abs(float(ebeam)) - abs(float(hms_p))
                    hms_th_rad = np.deg2rad(float(hms_th))
                    q2 = 4 * abs(float((ebeam)) * abs(float(hms_p)) * (np.sin(hms_th_rad/2))**2)
                    epsilon = 1 / (1 + 2 * (1 + (nu**2 / q2)) * np.tan(hms_th_rad/ 2))**2
                    xbj = q2 / (2 * mp * nu)

                    if USING_CURRENT_OFFSET:
                        current_offset = -0.0301
                        current_offset_corr = 1 / (1 + (current_offset / float(ibeam_2)))
                    else:
                        current_offset_corr = 1

                    weight = (float(boil_corr) * float(ps) * float(current_offset_corr)) / (float(livetime) * float(trackeff))

                    bigtable_rows_master.append({"runnum": runnum,
                                                 "run_type": run_type,
                                                 "target": target,
                                                 "phase": phase,
                                                 "xbj": xbj,
                                                 "q2": q2,
                                                 "ebeam": ebeam,
                                                 "ibeam_1": ibeam_1,
                                                 "ibeam_2": ibeam_2,
                                                 "ibeam_4a": ibeam_4a,
                                                 "ibeam_4c": ibeam_4c,
                                                 "qbeam_1": qbeam_1,
                                                 "qbeam_2": qbeam_2,
                                                 "qbeam_4a": qbeam_4a,
                                                 "qbeam_4c": qbeam_4c,
                                                 "hms_p": hms_p,
                                                 "hms_th": hms_th,
                                                 "hms_th_rad": hms_th_rad,
                                                 "ps3": ps3,
                                                 "ps4": ps4,
                                                 "trackeff": trackeff,
                                                 "nu": nu,
                                                 "boil_corr": boil_corr,
                                                 "epsilon": epsilon,
                                                 "weight": weight})
                except (ValueError, KeyError):
                    continue
bigtable_df_master = pd.DataFrame(bigtable_rows_master)

# -----------------------------------------------------
# Report file lookup
# -----------------------------------------------------
report_rows_master = []

for phase in phases:
    report_filepath = report_filepaths[phase]
    phase_runnums = [r["runnum"] for r in bigtable_rows_master if r["phase"] == phase]

    for runnum in phase_runnums:
        report_file = f"{report_filepath}/replay_hms_coin_production_{runnum}_-1.report"
        if not os.path.exists(report_file):
            continue
        threeoffour_rate = None
        threeoffour_counts = None
        elclean_rate = None
        elclean_counts = None

        with open(report_file) as rep_file:
            for repline in rep_file:
                if repline.startswith("HMS 3/4 Trigger Rate"):
                    m = re.search(r"\[\s*([0-9.]+)", repline)
                    if m:
                        threeoffour_rate = float(m.group(1))

                if repline.startswith("HMS 3/4 Trigger Rate"):
                    m = re.search(r":\s*([0-9.]+)", repline)
                    if m:
                        threeoffour_counts = float(m.group(1))

                if repline.startswith("hEL_CLEAN"):
                    m = re.search(r"\[\s*([0-9.]+)", repline)
                    if m:
                        elclean_rate = float(m.group(1))

                    n = re.search(r":\s*([0-9.]+)", repline)
                    if n:
                        elclean_counts = float(n.group(1))

            report_rows_master.append({"runnum": runnum,
                                       "3of4_rate": threeoffour_rate,
                                       "3of4_counts": threeoffour_counts,
                                       "elclean_rate": elclean_rate,
                                       "elclean_counts": elclean_counts})

report_df_master = pd.DataFrame(report_rows_master)
# -----------------------------------------------------
# Writing master file
# -----------------------------------------------------
master_df = bigtable_df_master.merge(report_df_master, on="runnum", how="left")

if not master_df.empty:
    master_df.to_csv(master_output_filepath, index=False)
    print(f"\n☢️  Saved {len(master_df)} total HMS DIS runs to {master_output_filepath}")
else:
    print(f"\n⚠️  No HMS DIS runs found for phase {phase}. No master file was written.")
            
# -----------------------------------------------------
# Filtering and writing out analysis tables
# -----------------------------------------------------
for phase in phases:
    if phase == "I":
        valid_beam_passes = ["4", "5"]
        keep_angle = {"4": vals["angle_4pass"],
                      "5": vals["angle_5pass"]}
    elif phase == "II":
        valid_beam_passes = ["3", "4", "5"]
        keep_angle = {"3": vals["angle_3pass_phaseII"],
                      "4": vals["angle_4pass_phaseII"],
                      "5": vals["angle_5pass_phaseII"]}

    for beam_pass in valid_beam_passes:
        _, beam_prefix = parse_beam_pass(beam_pass)
        for target in targets:
            filtered_df = master_df.copy()
            filtered_df = filtered_df[filtered_df["phase"] == phase]
            filtered_df = filtered_df[filtered_df["target"].str.lower() == target]
            filtered_df = filtered_df[filtered_df["ebeam"].astype(str).str.startswith(beam_prefix)]
            filtered_df = filtered_df[abs(filtered_df["hms_th"] - keep_angle[beam_pass]) <= angle_tolerance]

            setting = parse_setting(beam_pass, phase)
            output_dir = target.upper()
            output_filepath = f"{output_dir}/filtered_hmsdis_{setting}_{target}.csv"

            if not filtered_df.empty:
                filtered_df.to_csv(output_filepath, index=False)
                print(f"Saved {len(filtered_df)} runs to {output_filepath}")
