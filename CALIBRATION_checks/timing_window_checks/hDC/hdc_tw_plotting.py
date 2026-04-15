#!/usr/bin/env python3

import uproot
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import awkward as ak
from matplotlib.ticker import MultipleLocator, MaxNLocator

cuts = {"1u1": (-13300, -11275),
        "1u2": (-13300, -11250),
        "1x1": (-13300, -11225),
        "1x2": (-13300, -11175),
        "1v2": (-13400, -11275),
        "1v1": (-13200, -11025),
        "2v1": (-13300, -11150),
        "2v2": (-13200, -11000),
        "2x2": (-13300, -11025),
        "2x1": (-13300, -11125),
        "2u2": (-13400, -11175),
        "2u1": (-13200, -11125),}

pass0p1_min_cut = -13800.0
pass0p1_max_cut = -11400.0

if len(sys.argv) == 2:
    min_win = sys.argv[1]
else:
    min_win = input(f"Enter HDC TW min.  Options: 13<1,2,3,4,6,8>00, SCAN, FINAL:\n")

folder = f"/volatile/hallc/c-rsidis/relder/ROOTfiles/CALIB/HDC_TWMIN_{min_win}"
outdir = f"PNGs/{min_win}"
os.makedirs(outdir, exist_ok=True)

planes = ["1u1","1u2","2u1","2u2","1x1","1x2","2x1","2x2","1v1","1v2","2v1","2v2"]
files = sorted([f for f in os.listdir(folder) if f.endswith(".root")])
colors = plt.get_cmap("tab10").colors[:len(files)]

data_dict = {}
data_dict2 = {}
for fname in files:
    fullpath = os.path.join(folder, fname)
    with uproot.open(fullpath) as f:
        tree = f["T"]

        cer = tree["H.cer.npeSum"].array(library="ak")
        cal = tree["H.cal.etottracknorm"].array(library="ak")
        delta = tree["H.gtr.dp"].array(library="ak")

        cut = (cer > 3.0) & (cal > 0.9) & (np.abs(delta) < 8)

        data_dict[fname] = {}
        data_dict2[fname] = {}
        for p in planes:
            arr = tree[f"H.dc.{p}.dist"].array(library="ak")
            arr2 = tree[f"H.dc.{p}.rawtdc"].array(library="ak")
            data_dict[fname][p] = ak.flatten(arr[cut])
            data_dict2[fname][p] = ak.flatten(arr2[cut])
            # Without cuts, if needed...
            # data_dict[fname][p] = ak.flatten(arr)
            # data_dict2[fname][p] = ak.flatten(arr2[cut])
            
def short_name(fname):
    base = os.path.splitext(os.path.basename(fname))[0]
    parts = base.split("_")
    if len(parts) >= 3 and parts[0] == "hms" and parts[1] == "coin" and parts[2] == "replay" and parts[3] == "production":
        base = "_".join(parts[4:-2])
        return base
    return base

for plane in planes:
    fig, ax = plt.subplots(figsize=(14, 6))
    for i, fname in enumerate(files):
        data = ak.to_numpy(data_dict[fname][plane])
        ax.hist(data,bins=100,histtype="step",linewidth=1.5,label=short_name(fname),color=colors[i],)
    ax.set_xlabel(f"H.dc.{plane}.dist")
    ax.set_ylabel("Counts")
    ax.set_title(f"Time Window Testing, Plane {plane}")
    ax.legend(fontsize=8,loc="center left",bbox_to_anchor=(1.02, 0.5),frameon=True,)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"Plane_{plane}_min_{min_win}.png"), bbox_inches="tight")
    plt.close(fig)

for plane in planes:
    fig, ax = plt.subplots(figsize=(14, 6))
    for i, fname in enumerate(files):
        data = ak.to_numpy(data_dict2[fname][plane])
        ax.hist(data,bins=100,histtype="step",linewidth=1.5,label=short_name(fname),color=colors[i],)
    if plane in cuts:
        cut_min, cut_max = cuts[plane]
        ax.axvline(cut_min, color = "red", linestyle = "--", linewidth = 2, label = f"suggested cut min: {cut_min}, max: {cut_max}")
        ax.axvline(cut_max, color = "red", linestyle = "--", linewidth = 2)
    ax.axvline(pass0p1_min_cut, color = "black", linestyle = "--", linewidth = 2, label = f"pass0p1 cut min: {pass0p1_min_cut}, max: {pass0p1_max_cut}")
    ax.axvline(pass0p1_max_cut, color = "black", linestyle = "--", linewidth = 2)
    ax.set_xlabel(f"H.dc.{plane}.rawtdc")
    ax.set_ylabel("Counts")
    ax.xaxis.set_major_locator(MaxNLocator(nbins = 10))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    # ax.set_yscale("log")
    ax.set_title(f"Time Window Testing, Plane {plane}")
    ax.legend(fontsize=8,loc="center left",bbox_to_anchor=(1.02, 0.5),frameon=True,)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"Plane_{plane}_min_{min_win}_rawtdc.png"), bbox_inches="tight")
    plt.close(fig)



