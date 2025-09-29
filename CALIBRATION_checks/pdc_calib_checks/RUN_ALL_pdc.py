#!/usr/bin/env python3

import subprocess
from tqdm import tqdm

script_path = "CALIB_CHECKS_pdc.py"

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
        # print(result.stdout.decode())

    except subprocess.CalledProcessError as e:
        tqdm.write(f"Error running script for {runnum}: {e}")
        continue
    
