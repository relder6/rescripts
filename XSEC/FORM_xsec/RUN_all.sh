#!/usr/bin/bash

run_types=(hmsdis)
beam_passes=(4 5)
targets=(c cu ld2 lh2)

script_path="./FORM_xsec.py"

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

script_path_ratio="./FORM_ratios.py"
numerators=(c cu)
denominators=(ld2 lh2)

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for numerator in "${numerators[@]}"; do
            for denominator in "${denominators[@]}"; do
                echo
                echo "****************************************************"    
                echo "Running $run_type ratio at ${beam_pass}pass on ${numerator}/${denominator}"    
                "$script_path_ratio" "$run_type" "$beam_pass" "$numerator" "$denominator"
            done
        done
    done
done
