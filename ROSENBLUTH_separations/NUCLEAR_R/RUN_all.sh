#!/usr/bin/bash

run_types=(hmsdis)
beam_passes=(4)
targets=(c cu ld2 lh2)
nbins=1

script_path="./BIN_centering.py"

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for target in "${targets[@]}"; do
            echo
            echo "****************************************************"    
            echo "Running $run_type on $target"    
            "$script_path" "$run_type" "$beam_pass" "$target" "nbins"
        done
    done
done

run_types=(hmsdis)
beam_passes=(4)
targets=(c cu ld2 lh2)

script_path="./ROSENBLUTH_separation.py"

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for target in "${targets[@]}"; do
            echo
            echo "****************************************************"    
            echo "Running $run_type on $target"    
            "$script_path" "$run_type" "$beam_pass" "$target"
        done
    done
done
