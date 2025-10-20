#!/usr/bin/env python3

import subprocess

script_path = "ROOT_skim.py"

for runnum in range(23833, 25604):
    try:
        input_str = f"{runnum}"
        result = subprocess.run(
            ["python3", script_path],
            input = input_str.encode(),
            stdout = subprocess.DEVNULL,
            stderr = subprocess.DEVNULL,
            check = True
            )

    except subprocess.CalledProcessError as e:
        print(f"Error code {e.returncode} for run {runnum}; check if this is an HMS run.")
        continue
