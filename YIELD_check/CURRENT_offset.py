#!/usr/bin/env python3

import sys, os, csv
import matplotlib
import uproot
import numpy as np
import matplotlib.pyplot as plt
import boost_histogram as bh
import matplotlib.colors as mplcolors
from scipy.optimize import curve_fit

# Input file
csvfile = "CSVs/yield_check_hmsdis_4pass_c.csv"

data = np.genfromtxt(csvfile,delimiter=",",names=True,dtype=None,encoding=None)

polarity_mask = data["polarity"] == "-"

# runnum_mask = np.isin(data["runnum"], [24303, 24304, 24306, 24307, 24308])

runnum_mask = ~np.isin(data["runnum"], [0, 1])

targ_mask = np.isin(data["target"], ["cu", "c", "lh2"])

pass_mask = np.isin(data["beampass"], ["4Pass", "5Pass"])

current_mask = (data["current"] > 10) & (data["current"] < 45)

mask = polarity_mask & runnum_mask & current_mask & targ_mask & pass_mask

I = data["current"][mask]

yield_norm = data["yield"][mask]

yield_err = data["yield_err"][mask]

def fit(I, yield_corr, delta):
    return yield_corr * (1 + delta / I)

popt, pcov = curve_fit (fit, I, yield_norm, p0=[np.mean(yield_norm), 0.0] , sigma=yield_err, absolute_sigma=True)

yield_corr_fit, delta_fit = popt
yield_corr_err, delta_err = np.sqrt(np.diag(pcov))

print(f"Y_corr = {yield_corr_fit:.6f} ± {yield_corr_err:.6f}")
print(f"delta = {delta_fit:.6f} ± {delta_err:.6f}")

targs = np.unique(data["target"][mask])
passes = np.unique(data["beampass"][mask])
polarity = np.unique(data["polarity"][mask])

targ_str = ", ".join(targs)
pass_str = ", ".join(map(str, passes))
pol_str = ", ".join(polarity)

I_fit = np.linspace(min(I)*0.9, max(I)*1.1, 200)


yield_corrected = yield_norm / (1 + delta_fit / I)

y_top = (np.max(yield_norm) + np.max(yield_err))*1.1

fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize = (8.5, 11/2), gridspec_kw = {'height_ratios': [1,1]}, sharex = True, sharey = True)

ax_top.errorbar(I, yield_norm, yerr=yield_err, fmt='o', label='Unfitted Normalized Yields', color='red')
ax_top.plot(I_fit, fit(I_fit, yield_corr_fit, delta_fit), '-', label=rf'Fit: $\delta =$ {delta_fit:.4f}$\pm${delta_err:.4f} ', color='red')
ax_top.set_ylabel("Setting- and Charge-Normalized Yield")
ax_top.grid()
ax_top.legend()
ax_top.set_title(f"Current Offset / Luminosity Studies \nTargets: ({targ_str}), Pass: ({pass_str}), Polarity: ({pol_str})")

ax_bot.errorbar(I, yield_corrected, yerr=yield_err, fmt='o', label='Corrected Yield', color='black')
ax_bot.set_ylabel("Normalized Yield")
ax_bot.set_xlabel("Beam Current (µA)")
ax_bot.legend()
ax_bot.grid()

plt.tight_layout()
plt.show()
