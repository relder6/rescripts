This directory contains scripts to produce drift time : wire number plots for the HMS drift chamber.

It expects folders COIN/ and HMS/ within this directory, and XPLANE/ UPLANE/ and VPLANE/ directories beneath each.

It can be used for a single run,
[terminal] ./CALIB_CHECKS_hdc.py
[terminal] Input run number: <runnum goes here>

It can also be used to loop over many runs,
[terminal] ./RUN_ALL_hdc.py

Modify directory locations and shortnames for CALIB_CHECKS_hdc.py as needed.
Modify ranges for RUN_ALL_hdc.py as needed.

R. Elder <rmelder4929@gmail.com> Modified 30 Sept 2025
