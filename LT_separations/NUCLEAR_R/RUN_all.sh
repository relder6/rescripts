#!/usr/bin/bash

run_types=(hmsdis)
targets=(al c cu ld2 lh2)
phases=(i ii)
nbins=(1)

script_path="./BIN_centering.py"

for run_type in "${run_types[@]}"; do
    for target in "${targets[@]}"; do
        for phase in "${phases[@]}"; do
            
            echo
            echo "****************************************************"    
            echo "Bin Centering: $run_type on phase${phase} $target"    
            "$script_path" "$run_type" "$target" "${nbins[0]}" "$phase"
        done
    done
done

script_path="./ROSENBLUTH_separation.py"

for run_type in "${run_types[@]}"; do
    for target in "${targets[@]}"; do
        for phase in "${phases[@]}"; do
            echo
            echo "****************************************************"    
            echo "Performing Rosenbluth Separation of R: $run_type on phase${phase} $target"    
            "$script_path" "$run_type" "$target" "$phase"
        done
    done
done

