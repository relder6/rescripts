#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import boost_histogram as bh

target = 'Carbon'
run_number_electron = 24242
beam_charge = 18122.907  # uC
live_time = 1
tracking_eff = 1
prescale_factor = 1
run_number_positron = 24549
beam_charge_pos = 43022.260  # uC
live_time_pos = 1
tracking_eff_pos = 1
prescale_factor_pos = 3
normfac =  0.41089643900143341
USE_LOG_SCALE = False

beam_charge_milliC = beam_charge / 1000
beam_charge_milliC_pos = beam_charge_pos / 1000

# -----------------------------------------------------
# File paths
# -----------------------------------------------------

# directory_filepath = "/volatile/hallc/c-rsidis/cmorean/replay_pass0a/ROOTfiles"
# hms_pattern = f"hms_coin_replay_production_{runnum}_-1.root"

#electron_file = f"/w/hallc-scshelf2102/c-rsidis/relder/ROOTfiles/hms_coin_replay_production_{run_number_electron}_-1.root"
electron_file = f"/volatile/hallc/c-rsidis/cmorean/replay_pass0a/ROOTfiles/hms_coin_replay_production_{run_number_electron}_-1.root"
#positron_file = f"/w/hallc-scshelf2102/c-rsidis/relder/ROOTfiles/hms_coin_replay_production_{run_number_positron}_-1.root"
positron_file = f"/volatile/hallc/c-rsidis/cmorean/replay_pass0a/ROOTfiles/hms_coin_replay_production_{run_number_positron}_-1.root"
mc_file = f"/u/group/c-rsidis/relder/mc-single-arm/worksim/run_{run_number_electron}_mc_single_arm.root"

# -----------------------------------------------------
# Defining branches, using uproot to put them in data frames
# -----------------------------------------------------

branches = ["H.gtr.dp", "H.cal.etottracknorm", "H.gtr.ph", "H.gtr.th", "H.gtr.x", "H.gtr.y",
    "H.kin.Q2", "H.kin.x_bj", "H.kin.W", "H.cer.npeSum"]

branches_mc = ["hsxfp", "hsyfp", "hsxpfp", "hsypfp", "hsdelta", "hsytar", "hsxptar", "hsyptar",
    "hsztar", "fry", "xsnum", "ysnum", "xsieve", "ysieve", "stop_id", "xb", "q2",
    "w", "eprime", "theta", "sigvert", "sigrad", "weight"]

# -----------------------------------------------------
# Using uproot to load dataframes indexed over branches
# -----------------------------------------------------

df_elec = pd.DataFrame(uproot.open(electron_file)["T"].arrays(branches, library="np"))
print(f"Electron DF loaded: {len(df_elec)} rows")
df_pos = pd.DataFrame(uproot.open(positron_file)["T"].arrays(branches, library="np"))
print(f"Positron DF loaded: {len(df_pos)} rows")
#chunks = []
#for chunk in uproot.iterate(f"{mc_file}:h10", branches_mc, library="pd", step_size="100 MB"):
#    print(f"Loaded chunk with {len(chunk)} rows")
#    chunks.append(chunk)
#df_mc = pd.concat(chunks, ignore_index=True)
#print(f"MonteCarlo DF loaded: {len(df_mc)} rows")
#for chunk in uproot.iterate(f"{mc_file}:h10", branches_mc, library="pd", step_size="50 MB"):
#    print(f"Loaded chunk with {len(chunk)} rows")
#    df_mc = chunk
df_mc = pd.DataFrame(uproot.open(mc_file)["h10"].arrays(branches_mc, library="np"))
print(f"MonteCarlo DF loaded: {len(df_mc)} rows")

cut_elec = (df_elec["H.gtr.dp"].between(-8, 8) & (df_elec["H.cer.npeSum"] > 2) & (df_elec["H.cal.etottracknorm"] > 0.7))
cut_pos = (df_pos["H.gtr.dp"].between(-8, 8) & (df_pos["H.cer.npeSum"] > 2) & (df_pos["H.cal.etottracknorm"] > 0.7))
cut_mc = (df_mc["hsdelta"].between(-8, 8))

df_elec_cut = df_elec[cut_elec].copy()
df_pos_cut = df_pos[cut_pos].copy()
df_mc_cut = df_mc[cut_mc].copy()

# -----------------------------------------------------
# Weights
# -----------------------------------------------------

df_elec_cut["weights"] = prescale_factor / (beam_charge_milliC * tracking_eff * live_time)
df_pos_cut["weights"] = prescale_factor_pos / (beam_charge_milliC_pos * tracking_eff_pos * live_time_pos)
df_mc_cut["weights"] = df_mc_cut["weight"] * normfac * prescale_factor

# -----------------------------------------------------
# Binning
# -----------------------------------------------------

custom_bins = {
    "H.gtr.dp":            dict(num_bins=100, min=-10, max=10),
    "H.cal.etottracknorm": dict(num_bins=100, min=0.8, max=1.225),
    "H.gtr.ph":            dict(num_bins=200, min=-0.025, max=0.025),
    "H.gtr.th":            dict(num_bins=100, min=-0.075, max=0.075),
    "H.gtr.x":             dict(num_bins=100, min=-0.5, max=0.20),
    "H.gtr.y":             dict(num_bins=100, min=-1, max=1.25),
    "H.kin.Q2":            dict(num_bins=100, min=2.7, max=3.9),
    "H.kin.x_bj":          dict(num_bins=100, min=0.2, max=0.3),
    "H.kin.W":             dict(num_bins=100, min=3.15, max=3.4),
    "H.cer.npeSum":        dict(num_bins=100, min=0, max=20),
}

wide_bins = {k: dict(num_bins=100, min=-100, max=100) for k in custom_bins}

variable_mc_map = {
    "H.gtr.dp": "hsdelta",
    # "H.gtr.ph": "hsphi", <---- still working on the mapping!
    # "H.gtr.th": "hsth",
    # "H.gtr.x": "hsx",
    # "H.gtr.y": "hsy",
    "H.kin.Q2": "q2",
    "H.kin.x_bj": "xb",
    "H.kin.W": "w"
}

# --- Functions ---
def compute_ratio_with_uncertainty(num, den, var_num, var_den):
    ratio = np.full_like(num, np.nan, dtype=float)
    ratio_err = np.full_like(num, np.nan, dtype=float)
    valid = (den > 0) & (num > 0)
    ratio[valid] = num[valid] / den[valid]
    rel_err_num = np.zeros_like(num, dtype=float)
    rel_err_den = np.zeros_like(num, dtype=float)
    rel_err_num[valid] = np.sqrt(var_num[valid]) / num[valid]
    rel_err_den[valid] = np.sqrt(var_den[valid]) / den[valid]
    ratio_err[valid] = ratio[valid] * np.sqrt(rel_err_num[valid] ** 2 + rel_err_den[valid] ** 2)
    return ratio, ratio_err

def make_subtracted_vs_mc_figure(df_elec, df_pos, df_mc, variable, bin_info, use_log_scale=False):
    mc_var = variable_mc_map.get(variable)
    if mc_var is None or mc_var not in df_mc.columns:
        print(f"Warning: MC variable for '{variable}' not found or missing in MC DataFrame; skipping")
        return None, None, None, None, None, None, None, None

    axis = bh.axis.Regular(bin_info["num_bins"], bin_info["min"], bin_info["max"])
    bin_centers = axis.centers

    hist_elec = bh.Histogram(axis, storage=bh.storage.Weight())
    hist_pos = bh.Histogram(axis, storage=bh.storage.Weight())
    hist_mc = bh.Histogram(axis, storage=bh.storage.Weight())

    hist_elec.fill(df_elec_cut[variable].values, weight=df_elec_cut["weights"].values)
    hist_pos.fill(df_pos_cut[variable].values, weight=df_pos_cut["weights"].values)
    hist_mc.fill(df_mc_cut[mc_var].values, weight=df_mc_cut["weights"].values)

    vals_elec = hist_elec.view()['value']
    vars_elec = hist_elec.view()['variance']
    vals_pos = hist_pos.view()['value']
    vars_pos = hist_pos.view()['variance']
    vals_mc = hist_mc.view()['value']
    vars_mc = hist_mc.view()['variance']

    vals_sub = vals_elec - vals_pos
    vars_sub = vars_elec + vars_pos

    ratio, ratio_err = compute_ratio_with_uncertainty(vals_sub, vals_mc, vars_sub, vars_mc)

    total_sub = np.sum(vals_sub)
    total_sub_err = np.sqrt(np.sum(vars_sub))
    total_mc = np.sum(vals_mc)
    total_mc_err = np.sqrt(np.sum(vars_mc))

    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1, figsize=(8, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )

    ax_top.errorbar(
        bin_centers, vals_sub, yerr=np.sqrt(vars_sub), fmt="o", mfc="none",
        mec="green", ecolor="green", capsize=3, label="Electron - Positron (Data)"
    )
    ax_top.errorbar(
        bin_centers, vals_mc, yerr=np.sqrt(vars_mc), fmt="o", mfc="none",
        mec="purple", ecolor="purple", capsize=3, label="Monte Carlo"
    )
    ax_top.set_ylabel("Weighted Counts")
    ax_top.set_title(f"{variable} - Subtracted Yield vs MC & Ratio")
    ax_top.legend()
    ax_top.grid(True)
    if use_log_scale:
        ax_top.set_yscale("log")

    ax_bottom.errorbar(
        bin_centers, ratio, yerr=ratio_err, fmt="o", mfc="none",
        mec="black", ecolor="black", capsize=3, label="Subtracted / MC"
    )
    ax_bottom.axhline(1, color="gray", linestyle="--", linewidth=1)
    ax_bottom.set_xlabel(variable)
    ax_bottom.set_ylabel("Ratio")
    ax_bottom.grid(True)
    ax_bottom.legend()
    ax_bottom.set_xlim(bin_info["min"], bin_info["max"])

    plt.tight_layout()

    return fig, total_sub, total_sub_err, total_mc, total_mc_err, hist_elec, hist_pos, hist_mc

def save_all_subtracted_vs_mc_plots(pdf, df_elec, df_pos, df_mc, bins, use_log_scale=False):
    for variable, bin_info in bins.items():
        result = make_subtracted_vs_mc_figure(df_elec, df_pos, df_mc, variable, bin_info, use_log_scale)
        if result is None or result[0] is None:
            print(f"Skipping {variable}, no figure generated.")
            continue
        fig, *_ = result  # fixed unpacking
        pdf.savefig(fig)
        fig.savefig(f"{variable}_subtracted_vs_mc.png", dpi=300)
        plt.close(fig)

def compute_integrals_summary_subtracted_vs_mc(df_elec, df_pos, df_mc, bins):
    summary = []
    for variable, bin_info in bins.items():
        result = make_subtracted_vs_mc_figure(df_elec, df_pos, df_mc, variable, bin_info, use_log_scale=False)
        if result[0] is None:
            continue
        _, sub_tot, sub_err, mc_tot, mc_err, *_ = result
        if sub_tot > 0 and mc_tot > 0:
            ratio_tot = sub_tot / mc_tot
            ratio_tot_err = ratio_tot * np.sqrt((sub_err / sub_tot) ** 2 + (mc_err / mc_tot) ** 2)
        else:
            ratio_tot, ratio_tot_err = np.nan, np.nan
        summary.append({
            "Variable": variable,
            "Subtracted Total": sub_tot,
            "Subtracted Err": sub_err,
            "MC Total": mc_tot,
            "MC Err": mc_err,
            "Ratio of Totals": ratio_tot,
            "Ratio of Totals Err": ratio_tot_err
        })
    return summary

def add_summary_table_pdf_subtracted_vs_mc(pdf, integrals_summary, target, run_number_electron, run_number_positron):
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.axis("off")
    col_labels = ["Variable", "Subtracted Total ± Err", "MC Total ± Err", "Ratio of Totals ± Err"]
    table_data = [
        [
            entry["Variable"],
            f"{entry['Subtracted Total']:.4f} ± {entry['Subtracted Err']:.4f}",
            f"{entry['MC Total']:.4f} ± {entry['MC Err']:.4f}",
            f"{entry['Ratio of Totals']:.8f} ± {entry['Ratio of Totals Err']:.8f}"
            if not np.isnan(entry["Ratio of Totals"]) else "N/A",
        ]
        for entry in integrals_summary
    ]
    table = ax.table(cellText=table_data, colLabels=col_labels, loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)
    ax.set_title(
        f"Positron-Subtracted Yields vs MC for {target}, e- Run {run_number_electron} and e+ Run {run_number_positron}",
        fontsize=14, pad=20,
    )
    pdf.savefig(fig)
    plt.close(fig)

# --- Run ---
with PdfPages(f"Yield_Comparisons_{target}.pdf") as pdf:
    save_all_subtracted_vs_mc_plots(pdf, df_elec_cut, df_pos_cut, df_mc_cut, custom_bins, use_log_scale=USE_LOG_SCALE)
    integrals_summary = compute_integrals_summary_subtracted_vs_mc(df_elec_cut, df_pos_cut, df_mc_cut, wide_bins)
    add_summary_table_pdf_subtracted_vs_mc(pdf, integrals_summary, target, run_number_electron, run_number_positron)

print(f"Saved Yield_Comparisons_{target}.pdf with plots and integral summary table")
