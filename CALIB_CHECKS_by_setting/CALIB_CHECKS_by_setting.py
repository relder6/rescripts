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
from matplotlib.lines import Line2D
from tqdm import tqdm

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

input_settings_filepath = f"../FILTER_type/{selected_target_shortname.upper()}/{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_runs.dat" # The name and location of the lookup tables for the input runs by setting

runnums = []

with open(input_settings_filepath, "r") as infile:
    next(infile) # Skipping the header line here
    for line in infile:
        parts = line.strip().split("\t")
        runnums.append(parts[0])

root_directory = f"/volatile/hallc/c-rsidis/cmorean/replay_pass0a/ROOTfiles"
# coin_pattern = f"coin_replay_production_{runnum}_-1.root"
hms_pattern = f"hms_coin_replay_production_{runnum}_-1.root"

for runnum in runnums:
    try:
        input_root_filepath = f"{root_directory}/{hms_pattern}"
        f = uproot.open(input_root_filepath)
        file_type = "HMS"
        if not f or f.IsZombie():
            raise OSError(f"Failed to open file {input_root_filepath}")
    except OSError:
        print(f"Run {runnum} not valid for hdc calibration; skipping...")
        sys.exit(1)

print(f"Found {file_type} run {runnum}.  Processing...")

branches = ["H.gtr.dp", "H.cer.npeSum", "H.cal.etottracknorm"]

for runnum in runnums:
    try:
        with uproot.open(input_root_filepath) as file:
            if "T" not in file:
                print(f"No 'T' tree in run {runnum}, skipping...")
                continue

            tree = file["T"]

            missing_branches = [branch for branch in branches if branch not in tree.keys()]
            if missing_branches:
                for branch in missing_branches:
                    print(f"Branch {branch} not found in run {runnum}, skipping...")
                    continue

            if not missing_branches:

                data = tree.arrays(branches, library = "np")

                cut_cal_plots = (
                    (data["H.gtr.dp"].between(-8,8)) &
                    (data["H.cer.npeSum"] > 1.5)
                    )

                


# delta = "H.gtr.dp"

# hist = bh.Histogram(bh.axis.Regular(100, -10, 10), storage = bh.storage.Weight())

# results = []

# print(f"Starting analysis for {len(runnums)} runs...")

# for idx, (runnum, ps3, ps4, charge, eff, livetime) in enumerate(tqdm(zip(runnums, prescale3_factors, prescale4_factors, beam_charges, tracking_effs, livetimes), total=len(runnums))):
#     ps3_val = float(ps3)
#     ps4_val = float(ps4)

#     if ps3_val == -999 or ps4_val == -999:
#         print(f"WARNING: run {runnum} has issues with prescale values.  Fix the lookup table for this setting!")
#         continue
#     elif ps3_val != -1 and ps4_val != -1:
#         print(f"WARNING: run {runnum} has both prescales set (ps3={ps3_val}, ps4={ps4_val}), using ps3.")
#     elif ps3_val != -1 and ps4_val == -1:
#         ps = ps3_val
#     elif ps4_val != -1 and ps3_val == -1:
#         ps = ps4_val
#     else:
#         print(f"Run {runnum} has no valid prescale (both -1), skipping...")
#         continue

        
#     # input_root_filepath = f"~/rsidis-2025/hallc_replay_rsidis/ROOTfiles/hms_coin_replay_production_{runnum}_-1.root"
#     input_root_filepath = f"/volatile/hallc/c-rsidis/cmorean/replay_pass0a/ROOTfiles/hms_coin_replay_production_{runnum}_-1.root"
#     try:
#         with uproot.open(input_root_filepath) as file:
#             # print(f"Starting analysis for run {runnum}, continuing...")
#             if "T" not in file:
#                 print(f"No 'T' tree in run {runnum}, skipping...")
#                 continue

#             tree = file["T"]

#             if delta not in tree.keys():
#                 print(f"Branch {delta} not found in run {runnum}, skipping...")
#                 continue
            
#             data = tree.arrays(branches, library = "np")

#             cut = (
#                 (data["H.cer.npeSum"] > 2) &
#                 (data["H.cal.etottracknorm"] > 0.8) &
#                 (data["H.gtr.dp"] > -8.0) &
#                 (data["H.gtr.dp"] < 8.0)
#                 )

#             weight = 1000 * ps / ( float(livetime) * float(charge) * float(eff))
#             weights = np.full_like(data[delta], weight, dtype = float)

#             hist.reset()
#             hist.fill(data[delta][cut], weight=weights[cut])

#             delta_yield = hist.sum().value
#             delta_yield_err = hist.sum().variance**0.5

#             print(f"Run {runnum}: yield = {delta_yield:.6g}, error = ±{delta_yield_err:.6g}")

#             results.append({"runnum": runnum, "yield": delta_yield, "yield_err": delta_yield_err})

#             # print(f"Filled histogram for run {runnum}, moving on...")

#     except FileNotFoundError:
#         print(f"Missing file for run {runnum}, skipping...")
            
# df = pd.DataFrame(results)

# # #################

# # -- Categorizing runs based on central momentum polarity for plotting with different shapes
# momentum_values = np.array([float(m) for m in hmsmomentum[:len(df)]])
# pos_mask = momentum_values > 0
# elec_mask = momentum_values <= 0

# parsed_ps3 = np.array([float(ps3) for ps3 in prescale3_factors[:len(df)]])
# x_indices = np.arange(len(df))

# # -- Categorizing runs based on prescale for plotting with different colors
# ps3_mask = (parsed_ps3 != -999) & (parsed_ps3 != -1)
# ps4_mask = ~ps3_mask

# plt.figure(figsize=(12,6))

# # Ps3 electrons
# plt.scatter(x_indices[elec_mask & ps3_mask], df["yield"][elec_mask & ps3_mask],
#             color='red', s=20, marker='o', label='Ps3 HMS 3/4 elec')
# plt.errorbar(x_indices[elec_mask & ps3_mask], df["yield"][elec_mask & ps3_mask],
#              yerr=df["yield_err"][elec_mask & ps3_mask], fmt='none', ecolor='red')

# # Ps4 electrons
# plt.scatter(x_indices[elec_mask & ps4_mask], df["yield"][elec_mask & ps4_mask],
#             color='blue', s=20, marker='o', label='Ps4 HMS ELREAL elec')
# plt.errorbar(x_indices[elec_mask & ps4_mask], df["yield"][elec_mask & ps4_mask],
#              yerr=df["yield_err"][elec_mask & ps4_mask], fmt='none', ecolor='blue')

# # Ps3 positrons
# plt.scatter(x_indices[pos_mask & ps3_mask], df["yield"][pos_mask & ps3_mask],
#             color='red', s=20, marker='x', label='Ps3 HMS 3/4 pos')
# plt.errorbar(x_indices[pos_mask & ps3_mask], df["yield"][pos_mask & ps3_mask],
#              yerr=df["yield_err"][pos_mask & ps3_mask], fmt='none', ecolor='red')

# # Ps4 positrons
# plt.scatter(x_indices[pos_mask & ps4_mask], df["yield"][pos_mask & ps4_mask],
#             color='blue', s=20, marker='x', label='Ps4 HMS ELREAL pos')
# plt.errorbar(x_indices[pos_mask & ps4_mask], df["yield"][pos_mask & ps4_mask],
#              yerr=df["yield_err"][pos_mask & ps4_mask], fmt='none', ecolor='blue')

# plt.xticks(x_indices, df["runnum"], rotation=45, fontsize = 10)
# plt.xlabel("Run Number", fontsize = 12)
# plt.ylabel("Delta Yield per mC", fontsize = 12)
# plt.title(f"{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_yields", fontsize = 14)
# plt.grid(True, linestyle='--', color = 'gray', alpha=0.3)
# plt.margins(y = 0.1)

# plt.legend(frameon = True, fontsize = 10)

# plt.tight_layout()
# plt.savefig(f"{selected_type}_{selected_beam_pass}pass_{selected_target_shortname}_yields.png")


# ----- Below is copy of the script for hdc calibration check, will do similar looping process

#!/usr/bin/env python3

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kFatal
import sys

runnum = input(f"Input the run number you wish to analyze:")

# Directory information, modify as needed
root_directory = f"/volatile/hallc/c-rsidis/cmorean/replay_pass0a/ROOTfiles"
coin_pattern = f"coin_replay_production_{runnum}_-1.root"
hms_pattern = f"hms_coin_replay_production_{runnum}_-1.root"

try:
    input_root_filepath = f"{root_directory}/{coin_pattern}"
    f = ROOT.TFile.Open(input_root_filepath)
    file_type = "COIN"
    if not f or f.IsZombie():
        raise OSError(f"Failed to open file {input_root_filepath}")
except OSError:
    try:
        input_root_filepath = f"{root_directory}/{hms_pattern}"
        f = ROOT.TFile.Open(input_root_filepath)
        file_type = "HMS"
        if not f or f.IsZombie():
            raise OSError(f"Failed to open file {input_root_filepath}")
    except OSError:
        print(f"Run {runnum} not valid for hdc calibration; skipping...")
        sys.exit(1)
        
print(f"Found {file_type} run {runnum}.  Processing...")
        
def draw_group(tree, cards, cname, title):
    c = ROOT.TCanvas(cname, title, 1000, 800)
    c.Divide(2, 2, 0.01, 0.01)
    lines = []

    for i, cardnum in enumerate(cards):
        ROOT.gStyle.SetOptStat(1100)
        pad = c.cd(i+1)
        pad.SetLeftMargin(0.05)
        pad.SetRightMargin(0.05)
        pad.SetTopMargin(0.05)
        pad.SetBottomMargin(0.05)
        hname = f"h_{cardnum}"
        expr = f"H.dc.{cardnum}.time:H.dc.{cardnum}.wirenum >> {hname}(117,-5,110,110,-100,250)"
        tree.Draw(expr, "", "COLZ")
        h = ROOT.gPad.GetPrimitive(hname)

        stats = h.GetListOfFunctions().FindObject("stats")
        if stats:
            stats.SetTextSize(0.02)
            stats.SetX1NDC(0.75)
            stats.SetX2NDC(0.95)
            stats.SetY1NDC(0.75)
            stats.SetY2NDC(0.95)
            
        if h:
            # print(f"{hname} entries:", h.GetEntries()) # Uncomment this line if you wish to see number of entries printed for each card
            xmin = h.GetXaxis().GetXmin()
            xmax = h.GetXaxis().GetXmax()

            # Drawing a highlighted, dashed red line to appear at y = 0; this is to guide the eye while reviewing calibration plots.
            highlight = ROOT.TLine(xmin, 0, xmax, 0)
            highlight.SetLineColor(ROOT.kWhite)
            highlight.SetLineWidth(4)
            highlight.SetLineStyle(2)
            highlight.Draw("same")
            lines.append(highlight)
            
            line = ROOT.TLine(xmin, 0, xmax, 0)
            line.SetLineColor(ROOT.kRed)
            line.SetLineWidth(2)
            line.SetLineStyle(2)
            line.Draw("same")
            lines.append(line)
    c.lines = lines    
    c.Update()
    c.SaveAs(f"{cname}.png")
    return c

def main():
    T = f.Get("T")

    # Defining the cards
    us = ["1u1","1u2","2u1","2u2"]
    xs = ["1x1","1x2","2x1","2x2"]
    vs = ["1v1","1v2","2v1","2v2"]

    # Draw each card
    draw_group(T, xs, f"{file_type}/XPLANE/{file_type}_run_{runnum}_hdc_xplane", "X planes")
    draw_group(T, us, f"{file_type}/UPLANE/{file_type}_run_{runnum}_hdc_uplane", "U planes")
    draw_group(T, vs, f"{file_type}/VPLANE/{file_type}_run_{runnum}_hdc_vplane", "V planes")

if __name__ == "__main__":
    main()
