#!/usr/bin/env python3

import re

# Ask for input
input_run_number = input("Enter the run number: ")
input_filename = f"/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/REPORT_OUTPUT/HMS/PRODUCTION/output_get_good_dis_ev_{input_run_number}_-1.csv"
output_filename = "HMSDIS_stats.dat"

with open(input_filename, "r") as infile:
    lines = infile.readlines()

header_lines = []
stat_lines = []

for line in lines:
    if line.startswith("runnum"):
        header_lines.append(line.strip())
    else:
        stat_lines.append(line.strip())

# Extract and print header information
for header in header_lines:
    print(header)

# Process and print statistics
for stat in stat_lines:
    if stat:
        fields = re.split(r'\s*,\s*', stat)
        if len(fields) == 3:
            run_number = fields[0]
            good_dis = fields[1]
            norm_yield = fields[2]
            print(f"Run Number: {run_number}, Good Dis: {good_dis}, Norm Yield: {norm_yield}")
        else:
            print(f"Unexpected format in line: {stat}")

# Write output
with open(output_filename, "a") as outfile:
#    outfile.writelines(header_lines) # Don't need to rewrite this every row...
    outfile.writelines(stat_lines)
    outfile.write("\n")
