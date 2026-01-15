#!/usr/bin/env python3

import subprocess

run_types = ["HMSDIS"]
beam_passes = ["1", "4", "5"]
targets = ["al", "c", "cu", "ld2", "lh2", "dummy"]

script_path = "MAKE_csvs.py"

for run_type in run_types:
    for beam_pass in beam_passes:
        for target in targets:
            print(f"Calling script with run_type = {run_type}, beam_pass = {beam_pass}, target = {target}")

            input_str = f"{run_type}\n{beam_pass}\n{target}\n"

            result = subprocess.run(
                ["python3", script_path],
                input = input_str.encode(),
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE
                )

            print(result.stdout.decode())
            if result.stderr:
                print("Errors:", result.stderr.decode())
                
