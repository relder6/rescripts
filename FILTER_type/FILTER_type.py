#!/usr/bin/env python3

import os, re

input_filepath = "/w/hallc-scshelf2102/c-rsidis/relder/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat" # The location of the auxfiles runlist
# report_filepath = "/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/REPORT_OUTPUT/HMS/PRODUCTION/replay_hms_coin_production_{runnum}_-1.report"
report_filepath = "/work/hallc/c-rsidis/cmorean/replay_pass0a/REPORT_OUTPUT/HMS/PRODUCTION/replay_hms_coin_production_{runnum}_-1.report"


# -- User inputs and input processing--

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

# selected_angle = input("Enter the desired HMS angle (present options: ): ").strip() #  This is kind of pointless for HMSDIS, only one setting... commenting out

selected_target = input("Enter desired target (options: C, Cu, Al, LD2, LH2, Dummy): ").strip().lower()
selected_target_shortcut_to_target_variable = {"al":"aluminum","al13":"aluminum","aluminum":"aluminum",
                                               "c":"carbon","c12":"carbon","carbon":"carbon",
                                               "cu":"copper","cu29":"copper","copper":"copper",
                                               "opt1":"optics1","optics1":"optics1",
                                               "opt2":"optics2","optics2":"optics2",
                                               "d2":"ld2","ld2":"ld2",
                                               "h2":"lh2","lh2":"lh2",
                                               "hole":"hole","chole":"hole","c-hole":"hole",
                                               "dummy":"dummy","dum":"dummy",
                                               }
selected_target_shortname = selected_target_shortcut_to_target_variable.get(selected_target)
if not selected_target_shortname:
    print(f"Unknown target: {selected_target}.  Please try again.")
    exit(1)

output_filepath = f"{selected_target_shortname.upper()}/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_runs.dat" # The name and location of the output file

with open(input_filepath, "r") as infile:
    lines = infile.readlines()


# header_lines = [] # Not going to bother saving the header lines from the runlist.
run_lines = []

for line in lines:
    if line.lstrip().startswith("#"):        # dropping comment lines
        continue
    if line.lstrip().startswith("!"):        # dropping lines starting with !
        continue
    if re.match(r"^\s*[-=*]", line):         # dropping separator lines
        continue
    run_lines.append(line)

def extract_fields(line):
    parts = re.split(r'\s+', line.strip())
    ebeam = ""
    target_type = ""
    run_type = ""

    for part in parts:
        if re.match(r'^\d+\.\d+$', part): # looking for the first float in the line
            ebeam = part
            break

    run_types = ["hole", "optics", "heep", "hee", "hmsdis", "shmsdis", "pi-sidis", "pi+sidis", "junk"] # defining the known run types
    for part in parts:
        if part.lower() in run_types:
            run_type = part.lower()
            break

    targets = ["lh2", "carbon", "copper", "aluminum", "ld2", "dummy", "c-hole"] # defining the known targets
    for part in parts:
        if part.lower() in targets:
            target_type = part.lower()
            break

    return ebeam, target_type, run_type


filtered_lines = []

# debug_output_filepath = "debug_extract_fields.txt"
# debug_outfile = open(debug_output_filepath, "w")

# debug_output_filepath = "debug_extract_fields.txt"
# with open(debug_output_filepath, "w") as debug_outfile:
#     for line in run_lines:
#         ebeam, target_type, run_type = extract_fields(line)
#         debug_outfile.write(f"ebeam: {ebeam} target: {target_type} run_type: {run_type}\n")

for line in run_lines:
    ebeam, target_type, run_type = extract_fields(line)
    beam_match = ebeam.startswith(beam_prefix)
    if (
            selected_type == run_type.strip().lower() and
            selected_target_shortname == target_type.strip().lower() and
            beam_match
    ):
        filtered_lines.append(line)

# -- Adding function to make csv of runnums.  Then, reading it in, so runnums can be used below.
        
output_runnums_filepath = f"RUNNUMS/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_runnums.csv"

with open(output_runnums_filepath, "w") as outfile:
    runnums = [re.split(r'\s+',line.strip())[0] for line in filtered_lines]
    outfile.write(",".join(runnums))

with open(output_runnums_filepath, "r") as infile:
    runnums = infile.read().strip().split(",")

# -- Function to write the output file
        
with open(output_filepath, "w") as outfile:
    # Write header
    outfile.write("#Run#\tDate\ttStart\tEbeam\tIbeam\tTarget\tHMSp\tHMSth\tSHMSp\tSHMSth\tPrescaleSettings\tRunType\tBCM2CutCh\tPs3\tPs4\ttLive\tPTrigs\tELREAL\tTrackEf\tWeight\tnu\tQ2\tepsilon\t# Comments\n")
    
    for line in filtered_lines:
        parts = re.split(r'\s+', line.strip())
        # First 12 columns: standard columns
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
        comment = " ".join(parts[12:]) if len(parts) > 12 else "" #This joins everything after the 12th column together as a comment.
        if comment and not comment.startswith("# "):
            comment = "# " + comment
            
        # -- Extracting values from the report files

        report_file = report_filepath.format(runnum=runnum)
        bcm2cutch, ps3, ps4, trackeff, phystriggers, helreal = "-999", "-999", "-999", "-999", "-999", "-999"

        if os.path.exists(report_file):
            with open(report_file, "r") as rep:
                for rep_line in rep:
                    if rep_line.startswith("BCM2  Beam Cut Charge: "):
                        m = re.search(r"([\d.]+)\s+uC", rep_line)
                        if m:
                            bcm2cutch = m.group(1)
                    elif rep_line.startswith("Ps3_factor = "):
                        m = re.search(r"Ps3_factor\s*=\s*([+-]?\d+)", rep_line)
                        if m:
                            ps3 = float(m.group(1))
                    elif rep_line.startswith("Ps4_factor = "):
                        m = re.search(r"Ps4_factor\s*=\s*([+-]?\d+)", rep_line)
                        if m:
                            ps4 = float(m.group(1))
                    elif rep_line.startswith("E SING FID TRACK EFFIC         :     "):
                        m = re.search(r"([\d.]+)\s", rep_line)
                        if m:
                            trackeff = float(m.group(1))
                    elif rep_line.startswith("Accepted Physics Triggers      : "):
                        m = re.search(r"([\d.]+)\s", rep_line)
                        if m:
                            phystriggers = float(m.group(1))
                    elif rep_line.startswith("hEL_REAL  :	"):
                        m = re.search(r"([\d.]+)\s", rep_line)
                        if m:
                            helreal = float(m.group(1))

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
        livetime_unformatted = float(ps) * float(phystriggers) / float(helreal)
        if livetime_unformatted > 1:
            livetime_unformatted = 1
        livetime = float(livetime_unformatted)
        weight = float(ps) / ((float(livetime) * float(trackeff)))
        nu = float(abs(ebeam)) - float(abs(hms_p))
        Q2 = 4 * float(ebeam) * float(abs(hms_p)) * (np.sin(np.deg2rad(hms_th)/2))**2
        epsilon = 1 / (1 + 2 * (1 + (nu**2 / Q2)) * np.tan(np.deg2rad(hms_th) / 2)**2)
        # Composing the tsv line
        # tsv_line = "\t".join([runnum, date, tstart, ebeam, ibeam, target, hms_p, hms_th, shms_p, shms_th, prescales, runtype, bcm2cutch, ps3, ps4, livetime, phystriggers, helreal, trackeff, nu, Q2,  comment ])
        tsv_line = "\t".join(map(str, [runnum, date, tstart, ebeam, ibeam, target, hms_p, hms_th, shms_p, shms_th, prescales, runtype, bcm2cutch, f"{ps3:+}", f"{ps4:+}", f"{livetime:.3f}", f"{phystriggers:.8g}", f"{helreal:.8g}", f"{trackeff:.4f}", f"{weight:+.6f}", f"{nu:.3f}", f"{Q2:.3f}", f"{epsilon:.3f}", comment]))

        outfile.write(tsv_line + "\n")

print(f"\n☢️  Wrote {len(filtered_lines)} matching lines to {output_filepath}")
print(f"\n☢️  Wrote matching run numbers to {output_runnums_filepath}")
