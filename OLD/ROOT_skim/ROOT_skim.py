#!/usr/bin/env python3

import ROOT
import os
import re

ROOT.gErrorIgnoreLevel = ROOT.kError

runnum = input("Input desired run number.")

root_directory = f"/cache/hallc/c-rsidis/analysis/replays/pass0"
# coin_pattern = f"coin_replay_production_{runnum}_-1.root"
# shms_pattern = f"shms_coin_replay_production_{runnum}_-1.root"
hms_pattern = f"hms_coin_replay_production_{runnum}_-1.root"

input_filepath = f"{root_directory}/{hms_pattern}"

output_filepath = f"/work/hallc/c-rsidis/relder/pass0/reSKIM/skimmed_hms_coin_replay_production_{runnum}_-1.root"

selected_branches = [
    "H.gtr.x", "H.gtr.y", "H.gtr.dp", "H.gtr.p", "H.gtr.ph", "H.gtr.th", "H.gtr.beta", "H.gtr.index",
    "H.dc.x_fp", "H.dc.y_fp", "H.dc.xp_fp", "H.dc.yp_fp", "H.cer.npeSum", "H.cal.etottracknorm",
    "H.react.x", "H.react.y", "H.react.z", "H.kin.Q2", "H.kin.x_bj", "H.kin.W"
]

infile = ROOT.TFile(input_filepath, "READ")
intree = infile.Get("T")

intree.SetBranchStatus("*", 0) # This sets all branches to zero

for branches in selected_branches:
    intree.SetBranchStatus(branches, 1) # This maintains writing of all the selected branches

outfile = ROOT.TFile(output_filepath, "RECREATE")
skim_tree = intree.CloneTree(0)

for i in range(intree.GetEntries()):
    intree.GetEntry(i)
    skim_tree.Fill()

skim_tree.Write()
outfile.Close()
infile.Close()

print(f"Saved skimmed ROOT file to {output_filepath}")
