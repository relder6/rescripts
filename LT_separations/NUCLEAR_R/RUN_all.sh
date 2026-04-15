#!/usr/bin/bash

run_types=(hmsdis)
targets=(al c cu ld2 lh2)
nbins=(1)

script_path="./BIN_centering.py"

for run_type in "${run_types[@]}"; do
    for target in "${targets[@]}"; do
        echo
        echo "****************************************************"    
        echo "Bin Centering: $run_type on $target"    
        "$script_path" "$run_type" "$target" "${nbins[0]}"
    done
done

run_types=(hmsdis)
targets=(al c cu ld2 lh2)

script_path="./ROSENBLUTH_separation.py"

for run_type in "${run_types[@]}"; do
    for target in "${targets[@]}"; do
        echo
        echo "****************************************************"    
        echo "Performing Rosenbluth Separation of R: $run_type on $target"    
        "$script_path" "$run_type" "$target"
    done
done
