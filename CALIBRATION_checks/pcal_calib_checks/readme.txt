This directory contains scripts to assess the quality and stability of SHMS calorimeter calibrations.

It expects folders COIN/ and SHMS/ within this directory.  Make these first.

CALIB_CHECKS_pcal.py can be used for a single run,
[terminal] ./CALIB_CHECKS_pcal.py
[terminal] Input run number: <runnum goes here>


RUN_ALL_pcal.py can be used to loop over many runs,
[terminal] ./RUN_ALL_pcal.py
Running this script will produce FIT_pcal_results.dat.  This table is read in and processed further in the next script.

PLOT_pcal_fit.py can then be used to produce e/p vs run number plots,
[terminal] ./PLOT_pcal_fit.py

Modify directory locations and shortnames within the scripts as needed.
Modify ranges for RUN_ALL_pcal.py as needed.
Modify name of output .png from PLOT_... as needed; there is only one hard-coded name.

R. Elder <rmelder4929@gmail.com> Modified 13 October 2025
