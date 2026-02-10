#!/usr/bin/env python3

import os, sys, re
import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import uproot

hms_file = f"../../../hallc_replay_rsidis/PARAM/HMS/GEN/h_reftime_cut_coindaq.param"
shms_file = f"../../../hallc_replay_rsidis/PARAM/SHMS/GEN/p_reftime_cut.param"
coin_trig_file = f"../../../hallc_replay_rsidis/PARAM/TRIG/tcoin.param"
shms_trig_file = f"../../../hallc_replay_rsidis/PARAM/TRIG/tshms.param"
hms_trig_file = f"../../../hallc_replay_rsidis/PARAM/TRIG/thms.param"
rootfile_dir = f"/lustre24/expphy/volatile/hallc/c-rsidis/cmorean/replay_pass0_1/ROOTfiles"

if len(sys.argv) == 3:
    run_mode = sys.argv[1].strip().lower()
    runnum = sys.argv[2].strip()

notes = []

headers = []

# =====================================================================
# 'coin' mode
# =====================================================================
if run_mode == "coin":

    rootfile_location = f"{rootfile_dir}/coin_replay_production_{runnum}_-1.root"

    headers.append(f"****************************************\nFound COIN Reftimes\n****************************************\n")
    
    wanted_exact = {"t_coin_trig_adcrefcut",
                    "phodo_adcrefcut",
                    "pngcer_adcrefcut",
                    "phgcer_adcrefcut",
                    "paero_adcrefcut",
                    "pcal_adcrefcut",
                    "phodo_tdcrefcut",
                    "t_coin_trig_tdcrefcut",
                    "pdc_tdcrefcut",
                    "hhodo_adcrefcut",
                    "hcer_adcrefcut",
                    "hcal_adcrefcut",
                    "hhodo_tdcrefcut",
                    "hdc_tdcrefcut",
                    }

    rename_map = {"t_coin_trig_adcrefcut": "pFADC_TREF_ROC2",
                  "phodo_adcrefcut": "pFADC_TREF_ROC2",
                  "pngcer_adcrefcut": "pFADC_TREF_ROC2",
                  "phgcer_adcrefcut": "pFADC_TREF_ROC2",
                  "paero_adcrefcut": "pFADC_TREF_ROC2",
                  "pcal_adcrefcut": "pFADC_TREF_ROC2",
                  "phodo_tdcrefcut": "pT1",
                  "t_coin_trig_tdcrefcut": "pT2",
                  "pdc_tdcrefcut": "pDCREF(1,2,...,10)",
                  "hhodo_adcrefcut": "hFADC_TREF_ROC1",
                  "hcer_adcrefcut": "hFADC_TREF_ROC1",
                  "hcal_adcrefcut": "hFADC_TREF_ROC2",
                  "hhodo_tdcrefcut": "hT2",
                  "hdc_tdcrefcut": "hDCREF(1,2,...,5)",
                  }

    files_used = (coin_trig_file, hms_file, shms_file)

    notes.append("hT1 not used in COIN mode")

    print_order = ["pFADC_TREF_ROC2", "pT1", "pT2", "pDCREF(1,2,...,10)", "hFADC_TREF_ROC1", "hT1", "hT2", "hDCREF(1,2,...,5)"]    

# =====================================================================
# 'shms' mode
# =====================================================================
if run_mode == "shms":

    rootfile_location = f"{rootfile_dir}/shms_coin_replay_production_{runnum}_-1.root"

    headers.append(f"****************************************\nFound SHMS_coin Reftimes\n****************************************\n")
    wanted_exact = {"t_shms_trig_adcrefcut",
                    "phodo_adcrefcut",
                    "pngcer_adcrefcut",
                    "phgcer_adcrefcut",
                    "paero_adcrefcut",
                    "pcal_adcrefcut",
                    "phodo_tdcrefcut",
                    "t_shms_trig_tdcrefcut",
                    "pdc_tdcrefcut",
                    }

    rename_map = {"t_shms_trig_adcrefcut": "pFADC_TREF_ROC2",
                  "phodo_adcrefcut": "pFADC_TREF_ROC2",
                  "pngcer_adcrefcut": "pFADC_TREF_ROC2",
                  "phgcer_adcrefcut": "pFADC_TREF_ROC2",
                  "paero_adcrefcut": "pFADC_TREF_ROC2",
                  "pcal_adcrefcut": "pFADC_TREF_ROC2",
                  "phodo_tdcrefcut": "pT1",
                  "t_shms_trig_tdcrefcut": "pT2",
                  "pdc_tdcrefcut": "pDCREF(1,2,...,10)",
                  }

    files_used = (shms_trig_file, shms_file)

    print_order = ["pFADC_TREF_ROC2", "pT1", "pT2", "pDCREF(1,2,...,10)"]

# =====================================================================
# 'hms' mode
# =====================================================================
if run_mode == "hms":

    rootfile_location = f"{rootfile_dir}/hms_coin_replay_production_{runnum}_-1.root"

    headers.append(f"****************************************\nFound HMS_coin Reftimes\n****************************************\n")
    
    wanted_exact = {"t_hms_trig_adcrefcut",
                    "t_hms_trig_tdcrefcut",
                    "hhodo_adcrefcut",
                    "hcer_adcrefcut",
                    "hcal_adcrefcut",
                    "hhodo_tdcrefcut",
                    "hdc_tdcrefcut",
                    }

    rename_map = {"t_hms_trig_adcrefcut": "hFADC_TREF_ROC1",
                  "t_hms_trig_tdcrefcut": "hT1",
                  "hhodo_adcrefcut": "hFADC_TREF_ROC1",
                  "hcer_adcrefcut": "hFADC_TREF_ROC1",
                  "hcal_adcrefcut": "hFADC_TREF_ROC1",
                  "hhodo_tdcrefcut": "hT2",
                  "hdc_tdcrefcut": "hDCREF(1,2,...,5)",
                  }

    files_used = (hms_trig_file, hms_file)

    print_order = ["hFADC_TREF_ROC1", "hT1", "hT2", "hDCREF(1,2,...,5)"]    
    
# =====================================================================
# Printing Logic
# =====================================================================
values = {}

for filename in files_used:
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(";") or "=" not in line:
                continue

            var, val = line.split("=", 1)
            var = var.strip()
            valstr = val.strip()
            if var not in wanted_exact:
                continue
            try:
                val = float(valstr)
            except ValueError:
                continue

            key = rename_map[var]
            values.setdefault(key, []).append(val)

final = {}
for k, v in values.items():
    if len(set(v)) != 1:
        raise ValueError(f"Inconsistent values for {k}: {v}")
    final[k] = abs(v[0])

for header in headers:
    print(f"{header}")
    
for key in print_order:
    print(f"{key} = {final.get(key, '<not found>')}")

for note in notes:
    print(f"Note: {note}")

# =====================================================================
# Plotting Logic
# =====================================================================
pdf_name = f"PDFs/reftime_check_{run_mode.upper()}_run_{runnum}.pdf"

with uproot.open(rootfile_location) as f:
    tree = f["T"]

    with PdfPages(pdf_name) as pdf:
        simple_plots=[("pFADC_TREF_ROC2","T.coin.pFADC_TREF_ROC2_adcPulseTimeRaw","T.coin.pFADC_TREF_ROC2_adcMultiplicity"),
                      ("pT1","T.coin.pT1_tdcTimeRaw","T.coin.pT1_tdcMultiplicity"),
                      ("pT2","T.coin.pT2_tdcTimeRaw","T.coin.pT2_tdcMultiplicity"),
                      ("hFADC_TREF_ROC1","T.coin.hFADC_TREF_ROC1_adcPulseTimeRaw","T.coin.hFADC_TREF_ROC1_adcMultiplicity"),
                      ("hT2","T.coin.hT2_tdcTimeRaw","T.coin.hT2_tdcMultiplicity")]

        for key,leaf,multleaf in simple_plots:
            if key not in final: continue
            t=tree[leaf].array(library="np")
            m=tree[multleaf].array(library="np")
            m_nonzero=m[m>0]
            if len(m_nonzero)==0: continue
            # find the multiplicity that occurs most often
            vals,counts=np.unique(m_nonzero,return_counts=True)
            lowest_occ=np.min(m_nonzero)
            sel=(m==lowest_occ)&np.isfinite(t)
            fig,ax=plt.subplots(figsize=(8,5))
            ax.hist(t[sel],bins=200,histtype="step",linewidth=1.2)
            ax.axvline(final[key],color="red",linestyle="--",linewidth=1.5)
            ax.set_title(f"Run {runnum}: {leaf}\n{multleaf} = {lowest_occ}")
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("Counts")
            ax.grid(True,alpha=0.4)
            pdf.savefig(fig)
            plt.close(fig)

        def plot_dcref_all(prefix,nplanes,cut_value,title):
            fig,ax=plt.subplots(figsize=(8,5))
            best_mult_all=[]
            for i in range(1,nplanes+1):
                t=tree[f"T.coin.{prefix}{i}_tdcTimeRaw"].array(library="np")
                m=tree[f"T.coin.{prefix}{i}_tdcMultiplicity"].array(library="np")
                m_nonzero=m[m>0]
                if len(m_nonzero)==0: continue
                vals,counts=np.unique(m_nonzero,return_counts=True)
                lowest_occ=np.min(m_nonzero)
                sel=(m==lowest_occ)&np.isfinite(t)
                ax.hist(t[sel],bins=200,histtype="step",linewidth=1.1,label=f"{prefix}{i}")
                best_mult_all.extend(m[sel])
            if len(best_mult_all)==0: return
            median_mult=int(np.median(best_mult_all))
            ax.axvline(cut_value,color="red",linestyle="--",linewidth=1.5)
            dc_leaf = f"T.coin.{prefix}{{i}}_tdcTimeRaw"
            dc_mult_leaf = f"T.coin.{prefix}{{i}}_tdcMultiplicity"
            ax.set_title(f"Run {runnum}: {dc_leaf}\n{dc_mult_leaf} = {median_mult}")
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("Counts")
            ax.legend(fontsize=8,ncol=2)
            ax.grid(True,alpha=0.4)
            pdf.savefig(fig)
            plt.close(fig)

        plot_dcref_all("pDCREF",10,final["pDCREF(1,2,...,10)"],"SHMS DC Reference Time")
        plot_dcref_all("hDCREF",5,final["hDCREF(1,2,...,5)"],"HMS DC Reference Time")

print(f"Saved {pdf_name}.")
