TESTING_ROOTfiles - Developed by R. Elder, November 2025
!============================================

This folder contains two main scripts:
- EXTRACT_filenames.py: A script that will read a directory containing ROOTfiles, match them to desired patterns (i.e. 'shms_coin_replay_production_{number}_-1.root') and will write them to a comma separated list.  This list will be stored in filenames.txt located in the same directory in which you run this script.

- TESTING_ROOTfiles.py: A script that reads in the filenames.txt produced above and uses uproot to look for errors in certain branches.  You may alter the       branches as desired, but, more branches --> more time to run the script.  This script will produce error_log.txt, which contains an error report for all errors encountered.

Recommended use-case: verify transferred ROOTfiles
0) Transfer your ROOTfiles
1) Use EXTRACT_filenames.py in the original ROOTfile location to produce filenames.txt
2) Use TESTING_ROOTfiles.py to produce error_log.txt; consider re-naming.
3) Using the same filenames.txt, now use TESTING_ROOTfiles.py to produce an error_log.txt for the new directory.
4) Compare error logs.

