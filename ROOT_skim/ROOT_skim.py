#!/usr/bin/env python3

import ROOT
import os
import re

ROOT.gErrorIgnoreLevel = ROOT.kError

runnum = input("Input desired run number.")

root_directory = f"/volatile/hallc/c-rsidis/cmorean/replay_pass0a/ROOTfiles"
# coin_pattern = f"coin_replay_production_{runnum}_-1.root"
# shms_pattern = f"shms_coin_replay_production_{runnum}_-1.root"
hms_pattern = f"hms_coin_replay_production_{runnum}_-1.root"

input_filepath = f"{root_directory}/{hms_pattern}"

output_filepath = f"hms_skim_coin_replay_production_{runnum}_-1.root"

selected_branches = [
    "H.gtr.dp", "H.cal.etottracknorm", "H.gtr.ph", "H.gtr.th",
    "H.gtr.x", "H.gtr.y", "H.kin.Q2", "H.kin.x_bj",
    "H.kin.W", "H.cer.npeSum"
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
