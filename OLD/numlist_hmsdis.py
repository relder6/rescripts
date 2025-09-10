#!/usr/bin/env python3

import re

# Ask for input
input_filename = "/home/cdaq/rsidis-2025/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat"
selected_type = "HMSDIS"
output_filename = f"numlist_hmsdis.dat"

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
    for line in filtered_lines:
        fields = re.split(r'\s{2,}|\t+', line.strip())
        run_number = fields[0]
        outfile.write(run_number + ",")

print(f"\n✅ Wrote matching run numbers to {output_filename}")
