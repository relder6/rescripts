This directory contains scripts to produce drift time : wire number plots for the SHMS drift chamber.

It expects folders COIN/ and SHMS/ within this directory, and XPLANE/ UPLANE/ and VPLANE/ directories beneath each.

It can be used for a single run,
[terminal] ./CALIB_CHECKS_pdc.py
[terminal] Input run number: <runnum goes here>

It can also be used to loop over many runs,
[terminal] ./RUN_ALL_pdc.py

Modify directory locations and shortnames for CALIB_CHECKS_pdc.py as needed.
Modify ranges for RUN_ALL_pdc.py as needed.

R. Elder <rmelder4929@gmail.com> Modified 30 Sept 2025
