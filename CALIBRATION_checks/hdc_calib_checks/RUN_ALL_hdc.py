#!/usr/bin/env python3

import subprocess
from tqdm import tqdm

script_path = "CALIB_CHECKS_hdc.py"

for runnum in tqdm(range(23833, 25604)):
    try:
        input_str = f"{runnum}"
        result = subprocess.run(
            ["python3", script_path],
            input=input_str.encode(),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
            )
        # print(result.stdout.decode()) #Uncomment this line if you wish to see stdout prints.

    except subprocess.CalledProcessError as e:
        tqdm.write(f"Error code {e.returncode} for run {runnum}; check if this is an SHMS run.")
        continue
    
