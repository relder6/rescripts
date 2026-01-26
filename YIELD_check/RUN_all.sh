#!/usr/bin/bash

run_types=(hmsdis)
beam_passes=(4 5)
targets=(c cu ld2 lh2 dummy)

script_path="./YIELD_check.py"

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for target in "${targets[@]}"; do
            echo
            echo "****************************************************"
            echo "Running $run_type at ${beam_pass}pass on $target"
            "$script_path" "$run_type" "$beam_pass" "$target"
        done
    done
done
