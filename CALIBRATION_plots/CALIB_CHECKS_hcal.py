#!/usr/bin/env python3

import ROOT
ROOT.gROOT.SetBatch(True)

d_calo_fp = 292.64 #Distance, in centimeters, from the focal plane to the calorimeter position

xcalo = p.tr.x

def draw_group(tree, cards, cname, title):
    c = ROOT.TCanvas(cname, title, 1000, 800)
    c.Divide(2, 2)

    for i, cardnum in enumerate(cards):
        c.cd(i+1)
        hname = f"h_{cardnum}"
        expr = f"H.dc.{cardnum}.time:H.dc.{cardnum}.wirenum >> {hname}(122,-10,110,110,-100,350)"
        tree.Draw(expr, "", "COLZ")
        ROOT.gPad.SetRightMargin(0.15)
        h = ROOT.gPad.GetPrimitive(hname)
        if h:
            print(f"{hname} entries:", h.GetEntries())
    c.Update()
    c.SaveAs(f"{cname}.png")
    return c

def main():
    f = ROOT.TFile.Open("/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles/coin_replay_production_24240_-1.root")
    T = f.Get("T")  # adjust if tree has a different name

    # Define groups
    us = ["1u1","1u2","2u1","2u2"]
    xs = ["1x1","1x2","2x1","2x2"]
    vs = ["1v1","1v2","2v1","2v2"]

    # Draw each group
    draw_group(T, xs, "hdc_xplane", "X planes")
    draw_group(T, us, "hdc_uplane", "U planes")
    draw_group(T, vs, "hdc_vplane", "V planes")

if __name__ == "__main__":
    main()


