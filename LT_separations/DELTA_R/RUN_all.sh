#!/usr/bin/bash

run_types=(hmsdis)
num_targets=(al c cu)
denom_targets=(ld2)
nbins=(1)

script_path="./BIN_centering_model_ratios.py"

for run_type in "${run_types[@]}"; do
    for num_target in "${num_targets[@]}"; do
        for denom_target in "${denom_targets[@]}"; do
            echo
            echo "****************************************************"    
            echo "Bin Centering: $run_type on $num_target / $denom_target"    
            "$script_path" "$run_type" "$num_target" "$denom_target" "${nbins[0]}"
        done
    done
done

run_types=(hmsdis)
num_targets=(al c cu)
denom_targets=(ld2)

script_path="./ROSENBLUTH_separation.py"

for run_type in "${run_types[@]}"; do
    for num_target in "${num_targets[@]}"; do
        for denom_target in "${denom_targets[@]}"; do
            echo
            echo "****************************************************"    
            echo "Performing Rosenbluth Separation of delta R: $run_type on $num_target / $denom_target"    
            "$script_path" "$run_type" "$num_target" "$denom_target"
        done
    done
done
