#!/usr/bin/env python3

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kFatal
import sys

runnum = input(f"Input run number: ")

# Directory information, modify as needed
root_directory = f"/volatile/hallc/c-rsidis/cmorean/replay_pass0a/ROOTfiles"
coin_pattern = f"coin_replay_production_{runnum}_-1.root"
shms_pattern = f"shms_coin_replay_production_{runnum}_-1.root"

try:
    input_root_filepath = f"{root_directory}/{coin_pattern}"
    f = ROOT.TFile.Open(input_root_filepath)
    file_type = "COIN"
    if not f or f.IsZombie():
        raise OSError(f"Failed to open file {input_root_filepath}")
except OSError:
    try:
        input_root_filepath = f"{root_directory}/{shms_pattern}"
        f = ROOT.TFile.Open(input_root_filepath)
        file_type = "SHMS"
        if not f or f.IsZombie():
            raise OSError(f"Failed to open file {input_root_filepath}")
    except OSError:
        print(f"Run {runnum} not valid for pdc calibration; skipping...")
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
        expr = f"P.dc.{cardnum}.time:P.dc.{cardnum}.wirenum >> {hname}(117,-5,110,110,-100,250)"
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

    # Drawing each card
    draw_group(T, xs, f"{file_type}/XPLANE/{file_type}_run_{runnum}_pdc_xplane", "X planes")
    draw_group(T, us, f"{file_type}/UPLANE/{file_type}_run_{runnum}_pdc_uplane", "U planes")
    draw_group(T, vs, f"{file_type}/VPLANE/{file_type}_run_{runnum}_pdc_vplane", "V planes")

if __name__ == "__main__":
    main()


