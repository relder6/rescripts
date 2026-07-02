#!/usr/bin/bash

run_types=(hmsdis)
beam_passes=(1 3 4 5)
targets=(al c cu ld2 lh2 dummy)
phases=(i ii)

script_path="./FILTER_type.py"

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for target in "${targets[@]}"; do
            for phase in "${phases[@]}"; do
                "$script_path" "$run_type" "$beam_pass" "$target" "$phase"
            done
        done
    done
done


