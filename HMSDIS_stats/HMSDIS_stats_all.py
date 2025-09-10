#!/usr/bin/env python3

import re
import os
import csv

input_filename = "/home/cdaq/rsidis-2025/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat"
selected_type = "HMSDIS"
numlist_filename = "../TEMP/numlist_hmsdis.dat"
settings_filename = "../TEMP/settings_hmsdis.dat"
output_filename = "HMSDIS_stats_all.dat"

with open(input_filename, "r") as infile:
    lines = infile.readlines()

run_lines = []

for line in lines:
    if line.lstrip().startswith("#"):
        continue
    if line.lstrip().startswith("!"):
        continue
    if re.match(r"^\s*[-=*]", line):
        continue
    run_lines.append(line)

def extract_run_type(line):
    parts = re.split(r'\s{2,}|\t+', line.strip())
    if len(parts) >= 12:
        return parts[11]
    return ""

# Filtering the rsidis_runlist for only the selected HMSDIS type

filtered_lines = []

for line in run_lines:
    run_type = extract_run_type(line)
    if selected_type == run_type:
        filtered_lines.append(line)

# Write run numbers list, this will show up as a file in ../TEMP/
with open(numlist_filename, "w") as outfile:
    for line in filtered_lines:
        fields = re.split(r'\s{2,}|\t+', line.strip())
        run_number = fields[0]
        outfile.write(run_number + ",")
print(f"\n✅ Wrote matching run numbers to {numlist_filename}")

run_settings = {}
for line in filtered_lines:
    fields = re.split(r'\s{2,}|\t+', line.strip())
    run_number = fields[0]
    run_settings[run_number] = (
        fields[1] + "\t" + fields[2] + "\t" + fields[3] + "\t" +
        fields[4] + "\t" + fields[5] + "\t" + fields[6] + "\t" +
        fields[7] + "\t" + fields[8] + "\t" + fields[9] + "\t" +
        fields[10]
        )

with open(settings_filename, "w") as outfile:
    for line in filtered_lines:
        fields = re.split(r'\s{2,}|\t+', line.strip())
        run_number = fields[0]
        date = fields[1]
        start_time = fields[2]
        beam_energy = fields[3]
        target = fields[4]
        hms_p = fields[5]
        hms_th = fields[6]
        shms_p = fields[7]
        shms_th = fields[8]
        prescales = fields[9]
        comments = fields[10]
        outfile.write(
            run_number + "\t" + date + "\t" + start_time + "\t"
            + beam_energy + "\t" + target + "\t" + hms_p + "\t"
            + hms_th + "\t"
            #Commenting out the line below, don't need to have the shms details for this
            #+ shms_p + "\t" + shms_th + "\t"
            + prescales + "\t" + comments + "\n"
        )
print(f"\n✅ Wrote matching run numbers to {settings_filename}")

# Loop over run numbers to gather statistics
run_numbers = []
for line in filtered_lines:
    fields = re.split(r'\s{2,}|\t+', line.strip())
    run_number = fields[0]
    run_numbers.append(run_number)

# Removes the existing outfile if already in the directory
try:
    os.remove(output_filename)
except FileNotFoundError:
    print(f"No file found: {output_filename}, moving on...")
    
header_written = False
for run_number in run_numbers:
    stats_filename = f"/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/REPORT_OUTPUT/HMS/PRODUCTION/output_get_good_dis_ev_{run_number}_-1.csv"
    try:
        with open(stats_filename, "r") as infile:
            lines = infile.readlines()
            header_lines = []
            stat_lines = []
            for line in lines:
                if line.startswith("runnum"):
                    header_lines.append(line.strip())
                else:
                    stat_lines.append(line.strip())
            # Write header if not already done
            with open(output_filename, "a") as outfile:
                if not header_written and header_lines:
                    outfile.writelines([header + "\tdate\tstart_time\tbeam_energy\ttarget\thms_p\thms_th\tshms_p\tshms_th\tprescales\tcomments\n" for header in header_lines])
                    header_written = True
                # Write statistics lines
                for stat in stat_lines:
                    if stat:
                        settings_line = run_settings.get(run_number, "")
                        outfile.write(stat + "\t" + settings_line + "\n")
    except FileNotFoundError:
        print(f"Warning: statistics file missing for run {run_number}")

print(f"\n✅ Combined statistics written to {output_filename}")
