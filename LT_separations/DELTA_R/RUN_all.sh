#!/usr/bin/bash

run_types=(hmsdis)
num_targets=(al c cu)
denom_targets=(ld2)
phases=(i ii)
nbins=(1)

script_path="./BIN_centering_model_ratios.py"

for run_type in "${run_types[@]}"; do
    for num_target in "${num_targets[@]}"; do
        for denom_target in "${denom_targets[@]}"; do
            for phase in "${phases[@]}"; do
                
                echo
                echo "****************************************************"    
                echo "Bin Centering: $run_type on phase${phase} $num_target / $denom_target"    
                "$script_path" "$run_type" "$num_target" "$denom_target" "${nbins[0]}" "$phase"
            done
        done
    done
done

run_types=(hmsdis)
num_targets=(al c cu)
denom_targets=(ld2)
phases=(i ii)

script_path="./ROSENBLUTH_separation.py"

for run_type in "${run_types[@]}"; do
    for num_target in "${num_targets[@]}"; do
        for denom_target in "${denom_targets[@]}"; do
            for phase in "${phases[@]}"; do
                
                echo
                echo "****************************************************"    
                echo "Performing Rosenbluth Separation of delta R: $run_type on phase${phase} $num_target / $denom_target"    
                "$script_path" "$run_type" "$num_target" "$denom_target" "$phase"
            done
        done
    done
done

