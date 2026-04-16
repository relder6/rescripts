#!/usr/bin/env python3

import os, re
import numpy as np
import csv

# -----------------------------------------------------
# Defining paths here
# -----------------------------------------------------
input_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat" # The location of the auxfiles runlist
report_filepath = "/work/hallc/c-rsidis/replay/pass0p1/REPORT_OUTPUT/COIN/PRODUCTION/replay_coin_production_{runnum}_-1.report"
run_info_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_bigtable_pass0p1.csv"

output_filepath = "CSVs/coin_check.csv"

skip_runnums = []

# -----------------------------------------------------
# Reading auxfiles runlist, filtering and extracting lines
# -----------------------------------------------------
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

    # Only keep HMSDIS runs
    if (run_type == "pi+sidis" or run_type == "pi-sidis") and runnum not in skip_runnums:
        filtered_lines.append(line)

# -----------------------------------------------------
# Reading in big table values
# -----------------------------------------------------
runinfo_lookup = {}
if os.path.exists(run_info_filepath):
    with open(run_info_filepath, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                runnum = int(row["run"])
                runinfo_lookup[runnum] = {"bt_ebeam": row.get("ebeam", None),
                                          "bt_target": row.get("target", None),
                                          "bt_hms_p": row.get("hms_p", None),
                                          "bt_hms_th": row.get("hms_th", None),
                                          "bt_shms_p": row.get("shms_p", None),
                                          "bt_shms_th": row.get("shms_th", None),
                                          "bt_run_type": row.get("run_type", None),
                                          "bt_bcm2q": row.get("BCM2_Q", None),
                                          "bt_bcm2i": row.get("BCM2_I", None),
                                          "bt_ps1": row.get("ps1", None),
                                          "bt_ps2": row.get("ps2", None),
                                          "bt_ps3": row.get("ps3", None),
                                          "bt_ps4": row.get("ps4", None),
                                          "bt_ps5": row.get("ps5", None),
                                          "bt_ps6": row.get("ps6", None),
                                          "bt_comp_livetime": row.get("comp_livetime", None),
                                          "bt_boil_corr": row.get("boil_corr", None),
                                          "bt_fan_mean": row.get("fan_mean", None),
                                          "bt_Q2": row.get("Q2", None),
                                          "bt_xbj": row.get("x", None)}

            except (ValueError, KeyError) as e:
                print(f"Warning: skipping row due to error: {e}")

# -----------------------------------------------------
# Reading in report files and compiling the tsv
# -----------------------------------------------------
with open(output_filepath, "w") as outfile:
    outfile.write("runnum,ebeam_diff,ibeam_diff,targ,type,hms_p_diff,hms_th_diff,shms_p_diff,shms_th_diff,livetime_diff,ps5_diff,ps6_diff\n")

    for line in filtered_lines:
        parts = re.split(r'\s+', line.strip())
        runnum = parts[0]
        date = parts[1]
        tstart = parts[2]
        ebeam = parts[3]
        ibeam = parts[4]
        target = parts[5]
        hms_p = parts[6]
        hms_th = parts[7]
        shms_p = parts[8]
        shms_th = parts[9]
        prescales = parts[10]
        runtype = parts[11]
        comment = ""

        # Report file extraction
        report_file = report_filepath.format(runnum=runnum)
        h_bcm2cutch, p_bcm2cutch = "-999", "-999"
        ps5, ps6 = "-999", "-999"
        p_trackeff, h_trackeff = "-999", "-999"
        p_phystriggers, h_phystriggers = "-999", "-999"
        pelreal, helreal = "-999", "-999"

        if os.path.exists(report_file):
            with open(report_file, "r") as rep:
                for rep_line in rep:
                    
                    if rep_line.startswith("HMS BCM2  Beam Cut Charge: "):
                        m = re.search(r"([\d.]+)\s+uC", rep_line)
                        if m:
                            h_bcm2cutch = m.group(1)
                            
                    elif rep_line.startswith("SHMS BCM2  Beam Cut Charge: "):
                        m = re.search(r"([\d.]+)\s+uC", rep_line)
                        if m:
                            p_bcm2cutch = m.group(1)
                            
                    elif rep_line.startswith("Ps5_factor = "):
                        m = re.search(r"Ps5_factor\s*=\s*([+-]?\d+)", rep_line)
                        if m:
                            ps5 = float(m.group(1))
                            
                    elif rep_line.startswith("Ps6_factor = "):
                        m = re.search(r"Ps6_factor\s*=\s*([+-]?\d+)", rep_line)
                        if m:
                            ps6 = float(m.group(1))
                            
                    elif rep_line.startswith("E SING FID TRACK EFFIC"):
                        trackeff_count = 0
                        m = re.search(r"([\d.]+)\s", rep_line)
                        if m:
                            if trackeff_count == 0:
                                p_trackeff = float(m.group(1))
                            elif trackeff_count == 1:
                                h_trackeff = float(m.group(1))
                            trackeff_count += 1
                            
                    elif rep_line.startswith("SHMS Accepted Physics Triggers"):
                        m = re.search(r"([\d.]+)\s", rep_line)
                        if m:
                            p_phystriggers = float(m.group(1))
                            
                    elif rep_line.startswith("HMS Accepted Physics Triggers"):
                        m = re.search(r"([\d.]+)\s", rep_line)
                        if m:
                            h_phystriggers = float(m.group(1))
                            
                    elif rep_line.startswith("SHMS_pEL_REAL"):
                        m = re.search(r"([\d.]+)\s", rep_line)
                        if m:
                            pelreal = float(m.group(1))
                            
                    elif rep_line.startswith("HMS_hEL_REAL"):
                        m = re.search(r"([\d.]+)\s", rep_line)
                        if m:
                            helreal = float(m.group(1))

        ps5_val, ps6_val = float(ps5), float(ps6)
        if ps5_val == -999 or ps6_val == -999:
            print(f"REPORT FILE WARNING: Run {runnum} has issues with prescale values. Skipping...")
            continue
        elif ps5_val != -1 and ps6_val != -1:
            ps = ps5_val
        elif ps5_val != -1 and ps6_val == -1:
            ps = ps5_val
        elif ps6_val != -1 and ps5_val == -1:
            ps = ps6_val
        else:
            continue
        
        h_livetime_unformatted = float(ps) * float(h_phystriggers) / float(helreal)
        if h_livetime_unformatted > 1:
            h_livetime_unformatted = 1
        h_livetime = float(h_livetime_unformatted)

        p_livetime_unformatted = float(ps) * float(p_phystriggers) / float(pelreal)
        if p_livetime_unformatted > 1:
            p_livetime_unformatted = 1
        p_livetime = float(p_livetime_unformatted)
        
        # h_nu = abs(float(ebeam)) - abs(float(hms_p))
        # hms_th_rad = np.deg2rad(float(hms_th))
        # Q2 = 4 * abs(float(ebeam)) * abs(float(hms_p)) * (np.sin(hms_th_rad/2))**2
        # epsilon = 1 / (1 + 2 * (1 + (nu**2 / Q2)) * np.tan(hms_th_rad/2)**2)
        # xbj = Q2 / (2 * 0.938 * abs(float(nu)))

        runnum_int = int(runnum)
        runinfo_val = runinfo_lookup.get(runnum_int, {"bt_ebeam": None,
                                                      "bt_bcm2i": None,
                                                      "bt_run_type": None,
                                                      "bt_target": None,
                                                      "bt_fan_mean": None,
                                                      "bt_boil_corr": None,
                                                      "bt_hms_p": None,
                                                      "bt_hms_th": None,
                                                      "bt_shms_p": None,
                                                      "bt_shms_th": None,
                                                      "bt_comp_livetime": None,
                                                      "bt_ps5": None,
                                                      "bt_ps6": None,})

        # Bigtable extraction
        bt_ebeam = runinfo_val["bt_ebeam"]
        bt_bcm2i = runinfo_val["bt_bcm2i"]
        bt_run_type = runinfo_val["bt_run_type"]
        bt_target = runinfo_val["bt_target"]
        bt_Q2 = runinfo_val["bt_Q2"]
        bt_hms_p = runinfo_val["bt_hms_p"]
        bt_hms_th = runinfo_val["bt_hms_th"]
        bt_shms_p = runinfo_val["bt_shms_p"]
        bt_shms_th = runinfo_val["bt_shms_th"]
        bt_comp_livetime = runinfo_val["bt_comp_livetime"]
        bt_xbj = runinfo_val["bt_xbj"]
        bt_ps5 = runinfo_val["bt_ps5"]
        bt_ps6 = runinfo_val["bt_ps6"]

        try:
            ebeam_diff = float(ebeam) - float(bt_ebeam)
        except:
            ebeam_diff = None
        try:
            ibeam_diff = float(ibeam) - float(bt_bcm2i)
        except:
            ibeam_diff = None
        try:
            hms_p_diff = float(hms_p) - float(bt_hms_p)
        except:
            hms_p_diff = None
        try:
            hms_th_diff = float(hms_th) - float(bt_hms_th)
        except:
            hms_th_diff = None
        try:
            shms_p_diff = float(shms_p) - float(bt_shms_p)
        except:
            shms_p_diff = None
        try:
            shms_th_diff = float(shms_th) - float(bt_shms_th)
        except:
            shms_th_diff = None

        runtype_match = str(runtype).strip().lower() == str(bt_run_type).strip().lower()
        target_match = str(target).strip().lower() == str(bt_target).strip().lower()
        try:
            p_livetime_diff = float(p_livetime) - float(bt_comp_livetime)
        except:
            p_livetime_diff = None
        try:
            ps5_diff = float(ps5) - float(bt_ps5)
        except:
            ps5_diff = None
        try:
            ps6_diff = float(ps6) - float(bt_ps6)
        except:
            ps6_diff = None

        # Output CSV writing
        csv_line = ",".join(map(str, [runnum,
                                      f"{ebeam_diff:.6f}" if ebeam_diff is not None else "N/A",
                                      f"{ibeam_diff:.6f}" if ibeam_diff is not None else "N/A",
                                      target_match,
                                      runtype_match,
                                      f"{hms_p_diff:.6f}" if hms_p_diff is not None else "N/A",
                                      f"{hms_th_diff:.6f}" if hms_th_diff is not None else "N/A",
                                      f"{shms_p_diff:.6f}" if shms_p_diff is not None else "N/A",
                                      f"{shms_th_diff:.6f}" if shms_th_diff is not None else "N/A",
                                      f"{p_livetime_diff:.6f}" if p_livetime_diff is not None else "N/A",
                                      f"{ps5_diff:.6f}" if ps5_diff is not None else "N/A",
                                      f"{ps6_diff:.6f}" if ps6_diff is not None else "N/A",]))
        outfile.write(csv_line + "\n")

print(f"\n☢️  Wrote {len(filtered_lines)} HMSDIS runs to {output_filepath}")
