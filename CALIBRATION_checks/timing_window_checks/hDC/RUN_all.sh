#!/usr/bin/bash

inputs=(SCAN 13800 13600 13400 13300 13200 13100 13000 12900)

script_path="./hdc_tw_plotting.py"

for input in "${inputs[@]}"; do
    echo
    echo "***************************************"
    echo "Running hdc_tw min ${input}"
    "$script_path" "$input"
done
