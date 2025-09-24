#!/usr/bin/env python3

import re
import pandas as pd
import uproot
import boost_histogram as bh
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# -- User inputs and input processing, logic to select the setting--

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
    
input_runnums_filepath = f"../FILTER_type/RUNNUMS/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_runnums.csv" # The name and location of the input runnums file
input_settings_filepath = f"../FILTER_type/{selected_target_shortname.upper()}/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_runs.dat" # The name and location of the lookup tables for the input runs by setting


# with open(input_runnums_filepath, "r") as infile: # Realized I did NOT need to do this, I already have a convenient table with everything.  Duh!
#     runnums_line = infile.readline().strip()
#     runnums = [int(run) for run in runnums_line.split(",")]
#     print(f"Found list of run numbers for analysis: {runnums}")

# Leaving ordered list of all the components of the lookup table: runnum, date, tstart, ebeam, ibeam, target, hms_p, hms_th, shms_p, shms_th, prescales, runtype, bcm2cutch, ps3, ps4, livetime, comment

runnums = []
prescale3_factors = []
prescale4_factors = []
beam_charges = []
# tracking_effs = []
livetimes = []



with open(input_settings_filepath, "r") as infile:
    next(infile) # Skipping the header line here
    for line in infile:
        parts = line.strip().split("\t")
        runnums.append(parts[0])
        prescale3_factors.append(parts[13])
        prescale4_factors.append(parts[14])
        beam_charges.append(parts[12])
        livetimes.append(parts[15])

def parse_prescale(val):
    val = val.strip()
    if val in ("N/A", ""):
        return -999.0
    else:
        return float(val)

tracking_effs = [1.0] * len(runnums) # assuming tracking efficiencies are 1 for now

branches = ["H.gtr.dp", "H.cer.npeSum", "H.cal.etottracknorm"]

delta = "H.gtr.dp"

hist = bh.Histogram(bh.axis.Regular(100, -10, 10), storage = bh.storage.Weight())

results = []

for runnum, ps3, ps4, charge, eff, livetime in zip(
        runnums, prescale3_factors, prescale4_factors, beam_charges, tracking_effs, livetimes):
    ps3_val = parse_prescale(ps3)
    ps4_val = parse_prescale(ps4)

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

        
    input_root_filepath = f"~/rsidis-2025/hallc_replay_rsidis/ROOTfiles/hms_coin_replay_production_{runnum}_-1.root"
    try:
        with uproot.open(input_root_filepath) as file:
            
            if "T" not in file:
                print(f"No 'T' tree in run {runnum}, skipping...")
                continue

            tree = file["T"]

            if delta not in tree.keys():
                print(f"Branch {delta} not found in run {runnum}, skipping...")
                continue
            
            data = tree.arrays(branches, library = "np")

            cut = (
                (data["H.cer.npeSum"] > 1) &
                (data["H.cal.etottracknorm"] > 0.7)
                )

            weight = ps / ( float(livetime) * float(charge) * float(eff))
            weights = np.full_like(data[delta], weight, dtype = float)

            hist.reset()
            hist.fill(data[delta][cut], weight=weights[cut])

            delta_yield = hist.sum().value
            delta_yield_err = hist.sum().variance**0.5

            results.append({"runnum": runnum, "yield": delta_yield, "yield_err": delta_yield_err})

            print(f"Filled histogram for run {runnum}, moving on...")

    except FileNotFoundError:
        print(f"Missing file for run {runnum}, skipping...")
            
df = pd.DataFrame(results)

plt.figure(figsize=(12,6))

parsed_ps3 = [parse_prescale(ps3) for ps3 in prescale3_factors[:len(df)]]

colors = ['red' if ps_val != -999.0 and ps_val != -1 else 'blue' for ps_val in parsed_ps3]

plt.scatter(range(len(df)), df["yield"], c=colors, s=10)
plt.errorbar(range(len(df)), df["yield"], yerr=df["yield_err"], fmt='none', ecolor='black')

plt.xticks(range(len(df)), df["runnum"], rotation=45)
plt.xlabel("Run Number")
plt.ylabel("Yield")
plt.title(f"{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_yields")
plt.grid(True, linestyle='--', alpha=0.7)

red_patch = mpatches.Patch(color='red', label='Ps3 HMS 3/4')
blue_patch = mpatches.Patch(color='blue', label='Ps4 HMS ELREAL')
plt.legend(handles=[red_patch, blue_patch])

plt.tight_layout()
plt.savefig(f"{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_yields.png")
