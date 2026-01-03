#!/usr/bin/env python3

import re

# Ask for input
input_filename = "/home/cdaq/rsidis-2025/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat"
selected_type = "HMSDIS"
output_filename = f"runlist_filtered_by_hmsdis.dat"

with open(input_filename, "r") as infile:
    lines = infile.readlines()

run_lines = []

for line in lines:
    if line.lstrip().startswith("#"):        # drop comment lines
        continue
    if line.lstrip().startswith("!"):        # drop lines starting with !
        continue
    if re.match(r"^\s*[-=*]", line):         # drop separator lines
        continue
    run_lines.append(line)

# Pattern to match the Run_type field
def extract_run_type(line):
    parts = re.split(r'\s{2,}|\t+', line.strip())
    if len(parts) >= 12:
        return parts[11]
    return ""

filtered_lines = []

for line in run_lines:
    run_type = extract_run_type(line)
    if selected_type == run_type:
        filtered_lines.append(line)

# Write output
with open(output_filename, "w") as outfile:
    outfile.write("! Run#    Date       Start_Time   Ebeam I-beam Target   HMS_p   HMS_th  SHMS_p  SHMS_th PS1,2,3,4,5,6		Run_type	Comment\n")
    outfile.write("!---------------------------------------------------------------------------------------------------------------------------------------\n")
    outfile.writelines(filtered_lines)

print(f"\n✅ Wrote {len(filtered_lines)} matching lines to {output_filename}")
