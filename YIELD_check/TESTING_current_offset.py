#!/usr/bin/env python3

import sys, os, csv
import matplotlib
import uproot
import numpy as np
import matplotlib.pyplot as plt
import boost_histogram as bh
import matplotlib.colors as mplcolors
from scipy.optimize import curve_fit

selected_target = input(f"Input selected target, ld2 or lh2: ")

# Input file
csvfile = f"CSVs/yield_check_hmsdis_4pass_{selected_target}.csv"

data = np.genfromtxt(csvfile,delimiter=",",names=True,dtype=None,encoding=None)

polarity_mask = data["polarity"] == "-"

# runnum_mask = np.isin(data["runnum"], [24303, 24304, 24306, 24307, 24308, 24309])

runnum_mask = np.isin(data["runnum"], [24284,24285,24286,24287,24289,24290,24291,24292,24293,24294,24295,24296,24297,24298,24299])

targ_mask = np.isin(data["target"], ["cu", "c", "lh2", "ld2"])

pass_mask = np.isin(data["beampass"], ["4Pass", "5Pass"])

current_mask = (data["current"] > 5) & (data["current"] < 100)

mask = polarity_mask & runnum_mask & current_mask & targ_mask & pass_mask

current_offset = float(input("Input current offset to test: "))

I = data["current"][mask]

yield_norm = data["yield"][mask]

yield_err = data["yield_err"][mask]

yield_norm_corr = yield_norm / ( 1 + (current_offset / I))

yield_norm_corr_err = yield_err / (1 + (current_offset / I))

def fit(I, m, b):
    return m * I + b

popt, pcov = curve_fit(fit, I, yield_norm_corr, sigma=yield_norm_corr_err, absolute_sigma=True)

slope, intercept = popt
slope_err, int_err = np.sqrt(np.diag(pcov))

# yield_corr_fit, delta_fit = popt
# yield_corr_err, delta_err = np.sqrt(np.diag(pcov))

# print(f"Y_corr = {yield_corr_fit:.6f} ± {yield_corr_err:.6f}")
print(rf"yield_corr = ({slope:.4f} $\pm$ {slope_err:.4f}) * I + ({intercept:.4f} $\pm$ {int_err:.4f}) ")

I_fit = np.linspace(min(I)*0.9, max(I)*1.1, 200)

yield_change_per_100 = 100 * ((( slope * 100 + intercept ) - intercept) / (intercept))


# yield_corrected = yield_norm / (1 + delta_fit / I)

# y_top = (np.max(yield_norm) + np.max(yield_err))*1.1

fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize = (8.5, 11/2), gridspec_kw = {'height_ratios': [1,1]}, sharex = True, sharey = True)

ax_top.errorbar(I, yield_norm, yerr=yield_err, fmt='o', label='Unfitted Normalized Yields', color='navy')
ax_top.set_ylabel("Normalized Yield")
ax_top.grid()
ax_top.legend()
ax_top.set_title(rf"Current Offset Studies, Testing $\delta_I$ = {current_offset} on {selected_target.upper()}")

ax_bot.errorbar(I, yield_norm_corr, yerr=yield_norm_corr_err, fmt='o', label='Corrected Yield', color='orange')
ax_bot.plot(I_fit, fit(I_fit, slope, intercept), '--', label=rf"Fit: slope = {slope:.4f} $\pm$ {slope_err:.4f}""\n"rf"Yield Change per 100 $\mu$A = {yield_change_per_100:.4f}%", color='orange')
ax_bot.set_ylabel("Normalized Yield")
ax_bot.set_xlabel("Beam Current (µA)")
ax_bot.legend()
ax_bot.grid()

plt.tight_layout()
plt.show()
