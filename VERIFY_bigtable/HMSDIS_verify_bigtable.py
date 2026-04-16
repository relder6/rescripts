#!/usr/bin/env python3

import os, re
import numpy as np
import csv

# -----------------------------------------------------
# Defining paths here
# -----------------------------------------------------
input_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat"  # The location of the auxfiles runlist
report_filepath = "/work/hallc/c-rsidis/replay/pass0p1/REPORT_OUTPUT/HMS/PRODUCTION/replay_hms_coin_production_{runnum}_-1.report"
run_info_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_bigtable_pass0p1.csv"

output_filepath = "CSVs/hmsdis_check.csv"

skip_runnums = []


def safe_diff(a, b):
    try:
        return float(a) - float(b)
    except (TypeError, ValueError):
        return None

diff_warnings = {"ebeam_diff": 0.001,
                 "ibeam_diff": 0.001,
                 "qbeam_diff": 0.001,
                 "q2_diff": 0.001,
                 "xbj_diff": 0.001,
                 "hms_p_diff": 0.001,
                 "hms_th_diff": 0.0001,
                 "livetime_diff": 0.01,
                 "ps3_diff": 0.001,
                 "ps4_diff": 0.001}

warnings_by_run = {}

def check_warn(runnum, name, value):
    if value is None:
        return
    limit = diff_warnings.get(name)
    if limit is None:
        return
    if abs(value) > limit:
        entry = {"name": name, "value": float(value), "limit": float(limit)}
        warnings_by_run.setdefault(runnum, []).append(entry)
        
# -----------------------------------------------------
# Report file parsing
# -----------------------------------------------------
def parse_report(report_file):
    rf = {"rf_bcm2_cut_q": "-999",
          "rf_bcm2_cut_i": "-999",
          "rf_ps3": "-999",
          "rf_ps4": "-999",
          "rf_trackeff": "-999",
          "rf_phystriggers": "-999",
          "rf_helreal": "-999",
          "rf_3_of_4": "-999"}
    
    if not os.path.exists(report_file):
        return rf

    with open(report_file, "r") as rep:
        for rep_line in rep:
            if rep_line.startswith("BCM2  Beam Cut Charge: "):
                m = re.search(r"([\d.]+)\s+uC", rep_line)
                if m:
                    rf["rf_bcm2_cut_q"] = m.group(1)

            elif rep_line.startswith("BCM2 Beam Cut Current: "):
                m = re.search(r"([\d.]+)\s+uA", rep_line)
                if m:
                    rf["rf_bcm2_cut_i"] = m.group(1)

            elif rep_line.startswith("Ps3_factor = "):
                m = re.search(r"Ps3_factor\s*=\s*([+-]?\d+)", rep_line)
                if m:
                    rf["rf_ps3"] = float(m.group(1))

            elif rep_line.startswith("Ps4_factor = "):
                m = re.search(r"Ps4_factor\s*=\s*([+-]?\d+)", rep_line)
                if m:
                    rf["rf_ps4"] = float(m.group(1))

            elif rep_line.startswith("E SING FID TRACK EFFIC"):
                m = re.search(r"([\d.]+)\s", rep_line)
                if m:
                    rf["rf_trackeff"] = float(m.group(1))

            elif rep_line.startswith("Accepted Physics Triggers"):
                m = re.search(r"([\d.]+)\s", rep_line)
                if m:
                    rf["rf_phystriggers"] = float(m.group(1))

            elif rep_line.startswith("hEL_REAL"):
                m = re.search(r"([\d.]+)\s", rep_line)
                if m:
                    rf["rf_helreal"] = float(m.group(1))

            elif rep_line.startswith("HMS 3/4 Triggers"):
                  m = re.search(r":\s*(\d+)", rep_line)
                  if m:
                      rf["rf_3_of_4"] = float(m.group(1))

    return rf


# -----------------------------------------------------
# Auxfile parsing
# -----------------------------------------------------
with open(input_filepath, "r") as infile:
    lines = infile.readlines()

run_lines = []

for line in lines:
    if line.lstrip().startswith("#"): # dropping comment lines
        continue
    if line.lstrip().startswith("!"): # dropping lines starting with !
        continue
    if re.match(r"^\s*[-=*]", line): # dropping separator lines
        continue
    run_lines.append(line)

runlist_rows = []

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

    if (run_type == "hmsdis") and runnum not in skip_runnums:
        runlist_rows.append({"runnum": runnum,
                             "rl_date": parts[1],
                             "rl_tstart": parts[2],
                             "rl_ebeam": parts[3],
                             "rl_target": parts[5],
                             "rl_hms_p": parts[6],
                             "rl_hms_th": parts[7],
                             "rl_shms_p": parts[8],
                             "rl_shms_th": parts[9],
                             "rl_prescales": parts[10],
                             "rl_runtype": parts[11]})

# -----------------------------------------------------
# Bigtable parsing
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
                                          "bt_comp_livetime": row.get("comp_livetime", None),
                                          "bt_boil_corr": row.get("boil_corr", None),
                                          "bt_fan_mean": row.get("fan_mean", None),
                                          "bt_Q2": row.get("Q2", None),
                                          "bt_xbj": row.get("x", None),}
            except (ValueError, KeyError) as e:
                print(f"Warning: skipping row due to error: {e}")

# -----------------------------------------------------
# Now actually reading in all the values, computing differences, writing output file
# -----------------------------------------------------
with open(output_filepath, "w") as outfile:
    outfile.write("runnum,ebeam_diff,ibeam_diff,qbeam_diff,targ,type,q2_difference,xbj_difference,hms_p_diff,hms_th_diff,shms_p_diff,shms_th_diff,livetime_diff,ps3_diff,ps4_diff\n")

    for row in runlist_rows:
        runnum = row["runnum"]
        rl_ebeam = row["rl_ebeam"]
        rl_target = row["rl_target"]
        rl_hms_p = row["rl_hms_p"]
        rl_hms_th = row["rl_hms_th"]
        rl_shms_p = row["rl_shms_p"]
        rl_shms_th = row["rl_shms_th"]
        rl_runtype = row["rl_runtype"]

        report_file = report_filepath.format(runnum=runnum)
        rf = parse_report(report_file)

        ps3_val = float(rf["rf_ps3"])
        ps4_val = float(rf["rf_ps4"])
        if ps3_val == -999 or ps4_val == -999:
            print(f"REPORT FILE WARNING: Run {runnum} has issues with prescale values. Skipping...")
            continue
        if ps3_val != -1 and ps4_val != -1:
            ps = ps3_val
            livetime_unformatted = float(ps) * float(rf["rf_phystriggers"]) / float(rf["rf_3_of_4"])
            if livetime_unformatted > 1:
                livetime_unformatted = 1
            livetime = float(livetime_unformatted)
        elif ps3_val != -1 and ps4_val == -1:
            ps = ps3_val
            livetime_unformatted = float(ps) * float(rf["rf_phystriggers"]) / float(rf["rf_3_of_4"])
            if livetime_unformatted > 1:
                livetime_unformatted = 1
            livetime = float(livetime_unformatted)
        elif ps4_val != -1 and ps3_val == -1:
            ps = ps4_val
            livetime_unformatted = float(ps) * float(rf["rf_phystriggers"]) / float(rf["rf_helreal"])
            if livetime_unformatted > 1:
                livetime_unformatted = 1
            livetime = float(livetime_unformatted)
        else:
            continue



        ebeam = float(rl_ebeam)
        hms_p = float(rl_hms_p)
        hms_th = float(rl_hms_th)
        nu = abs(ebeam) - abs(hms_p)
        hms_th_rad = np.deg2rad(hms_th)
        Q2 = 4 * abs(ebeam) * abs(hms_p) * (np.sin(hms_th_rad / 2))**2
        epsilon = 1 / (1 + 2 * (1 + (nu**2 / Q2)) * np.tan(hms_th_rad / 2)**2)
        xbj = Q2 / (2 * 0.938 * abs(float(nu)))

        runnum_int = int(runnum)
        runinfo_val = runinfo_lookup.get(runnum_int, {})
        bt = {"bt_ebeam": runinfo_val.get("bt_ebeam"),
              "bt_bcm2i": runinfo_val.get("bt_bcm2i"),
              "bt_bcm2q": runinfo_val.get("bt_bcm2q"),
              "bt_run_type": runinfo_val.get("bt_run_type"),
              "bt_target": runinfo_val.get("bt_target"),
              "bt_Q2": runinfo_val.get("bt_Q2"),
              "bt_hms_p": runinfo_val.get("bt_hms_p"),
              "bt_hms_th": runinfo_val.get("bt_hms_th"),
              "bt_shms_p": runinfo_val.get("bt_shms_p"),
              "bt_shms_th": runinfo_val.get("bt_shms_th"),
              "bt_comp_livetime": runinfo_val.get("bt_comp_livetime"),
              "bt_xbj": runinfo_val.get("bt_xbj"),
              "bt_ps3": runinfo_val.get("bt_ps3"),
              "bt_ps4": runinfo_val.get("bt_ps4")}

        runtype_match = str(rl_runtype).strip().lower() == str(bt["bt_run_type"]).strip().lower()
        target_match = str(rl_target).strip().lower() == str(bt["bt_target"]).strip().lower()

        ebeam_diff = safe_diff(rl_ebeam, bt["bt_ebeam"])
        ibeam_diff = safe_diff(rf["rf_bcm2_cut_i"], bt["bt_bcm2i"])
        qbeam_diff = safe_diff(rf["rf_bcm2_cut_q"], bt["bt_bcm2q"])
        hms_p_diff = safe_diff(rl_hms_p, bt["bt_hms_p"])
        hms_th_diff = safe_diff(rl_hms_th, bt["bt_hms_th"])
        shms_p_diff = safe_diff(rl_shms_p, bt["bt_shms_p"])
        shms_th_diff = safe_diff(rl_shms_th, bt["bt_shms_th"])
        livetime_diff = safe_diff(livetime, bt["bt_comp_livetime"])
        Q2_diff = safe_diff(Q2, bt["bt_Q2"])
        xbj_diff = safe_diff(xbj, bt["bt_xbj"])
        ps3_diff = safe_diff(rf["rf_ps3"], bt["bt_ps3"])
        ps4_diff = safe_diff(rf["rf_ps4"], bt["bt_ps4"])

        check_warn(runnum, "ebeam_diff", ebeam_diff)
        check_warn(runnum, "ibeam_diff", ibeam_diff)
        check_warn(runnum, "qbeam_diff", qbeam_diff)
        check_warn(runnum, "hms_p_diff", hms_p_diff)
        check_warn(runnum, "hms_th_diff", hms_th_diff)
        check_warn(runnum, "shms_p_diff", shms_p_diff)
        check_warn(runnum, "shms_th_diff", shms_th_diff)
        check_warn(runnum, "livetime_diff", livetime_diff)
        check_warn(runnum, "Q2_diff", Q2_diff)
        check_warn(runnum, "xbj_diff", xbj_diff)
        check_warn(runnum, "ps3_diff", ps3_diff)
        check_warn(runnum, "ps4_diff", ps4_diff)

        csv_line = ",".join(map(str, [runnum,
                                      f"{ebeam_diff:.6f}" if ebeam_diff is not None else "N/A",
                                      f"{ibeam_diff:.6f}" if ibeam_diff is not None else "N/A",
                                      f"{qbeam_diff:.6f}" if qbeam_diff is not None else "N/A",
                                      target_match,
                                      runtype_match,
                                      f"{Q2_diff:.6f}" if Q2_diff is not None else "N/A",
                                      f"{xbj_diff:.6f}" if xbj_diff is not None else "N/A",
                                      f"{hms_p_diff:.6f}" if hms_p_diff is not None else "N/A",
                                      f"{hms_th_diff:.6f}" if hms_th_diff is not None else "N/A",
                                      f"{shms_p_diff:.6f}" if shms_p_diff is not None else "N/A",
                                      f"{shms_th_diff:.6f}" if shms_th_diff is not None else "N/A",
                                      f"{livetime_diff:.6f}" if livetime_diff is not None else "N/A",
                                      f"{ps3_diff:.6f}" if ps3_diff is not None else "N/A",
                                      f"{ps4_diff:.6f}" if ps4_diff is not None else "N/A"]))
        outfile.write(csv_line + "\n")
#print("-----------------------------------------------------")
print(f"SUMMARY OF FOUND WARNINGS: field, difference, limit")
#print("-----------------------------------------------------")
for runnum, entries in warnings_by_run.items():
    print(f"Run {runnum}--------------------------------------------")
    header_field = "field"
    header_diff = "diff"
    header_limit = "limit"

    w_field = 13
    w_diff = 13
    w_limit = 13

    #print(f"  {header_field:<{w_field}} {header_diff:>{w_diff}} {header_limit:>{w_limit}}")
    # print(f"  {'-'*w_field} {'-'*w_diff} {'-'*w_limit}")

    for e in entries:
        name = e["name"]
        val = e["value"]
        limit = e["limit"]
        print(f"  {name:<{w_field}} {val:>{w_diff}.6f} {limit:>{w_limit}.6f}")
    # print()

print(f"\n☢️  Wrote {len(runlist_rows)} HMSDIS runs to {output_filepath}")
