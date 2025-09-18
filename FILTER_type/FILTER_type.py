#!/usr/bin/env python3

import re

input_filename = "/home/cdaq/rsidis-2025/hallc_replay_rsidis/AUX_FILES/rsidis_runlist.dat" # The location of the file to be read in


# --User inputs and input processing--

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

# selected_angle = input("Enter the desired HMS angle (present options: ): ").strip() #  This isn't a good idea for HMSDIS...

selected_target = input("Enter desired target (options: C, Cu, Al, LD2, Dummy): ").strip().lower()
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

output_filename = f"{selected_target_shortname.upper()}/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_runs.dat" # The name and location of the output file

with open(input_filename, "r") as infile:
    lines = infile.readlines()

# Header lines
header_lines = []
run_lines = []

for line in lines:
    if line.lstrip().startswith("#"):        # drop comment lines
        continue
    if line.lstrip().startswith("!"):        # drop lines starting with !
        continue
    if re.match(r"^\s*[-=*]", line):         # drop separator lines
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

# debug_output_filename = "debug_extract_fields.txt"
# debug_outfile = open(debug_output_filename, "w")

# debug_output_filename = "debug_extract_fields.txt"
# with open(debug_output_filename, "w") as debug_outfile:
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

#expected_columns = 13  # Run# through Run_type (12) + Comment (13th)
with open(output_filename, "w") as outfile:
    # Write header
    outfile.write("#Run#\tdate\ttstart\tebeam\tIbeam\ttarget\tHMSp\tHMSth\tSHMSp\tSHMSth\tps1,ps2,ps3,ps4,ps5,ps6\truntype\t# comments\n")
    for line in filtered_lines:
        parts = re.split(r'\s+', line.strip())
        # First 12 columns: standard columns
        fixed = parts[:12]
        # Everything after the 12th part gets joined together as the comment
        comment = " ".join(parts[12:]) if len(parts) > 12 else ""
        comment_field = f"# {comment}" if comment else ""   # Only add marker if there is a comment
        # Compose the TSV line
        tsv_line = "\t".join(fixed + [comment_field])
        outfile.write(tsv_line + "\n")


#####


print(f"\n✅ Wrote {len(filtered_lines)} matching lines to {output_filename}")
