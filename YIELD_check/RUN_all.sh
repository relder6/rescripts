#!/usr/bin/bash

run_types=(hmsdis)
beam_passes=(3 4 5)
targets=(al c cu ld2 lh2)
phases=(i ii)

script_path="./YIELD_check.py"

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for target in "${targets[@]}"; do
            for phase in "${phases[@]}"; do
                echo
                echo "****************************************************"
                echo "Running $run_type at ${beam_pass}pass on $target phase $phase"
                "$script_path" "$run_type" "$beam_pass" "$target" "$phase"
            done
        done
    done
done

