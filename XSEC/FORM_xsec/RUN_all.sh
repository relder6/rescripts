#!/usr/bin/bash

run_types=(hmsdis)
beam_passes=(3 4 5)
targets=(al c cu ld2 lh2)
phases=(i ii)

script_path="./FORM_xsec.py"

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for target in "${targets[@]}"; do
            for phase in "${phases[@]}"; do

                if [[ "$beam_pass" == "3" && "$phase" == "i" ]]; then
                    continue
                fi
                
                echo
                echo "****************************************************"    
                echo "Running $run_type at ${beam_pass}pass on $target phase${phase}"    
                "$script_path" "$run_type" "$beam_pass" "$target" "$phase"
            done
        done
    done
done

script_path_ratio="./FORM_ratios.py"
numerators=(al c cu)
denominators=(ld2 lh2)
phases=(i ii)

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for numerator in "${numerators[@]}"; do
            for denominator in "${denominators[@]}"; do
                for phase in "${phases[@]}"; do

                    if [[ "$beam_pass" == "3" && "$phase" == "i" ]]; then
                        continue
                    fi
                    
                    echo
                    echo "****************************************************"    
                    echo "Running $run_type ratio at ${beam_pass}pass phase${phase} on ${numerator}/${denominator}"    
                    "$script_path_ratio" "$run_type" "$beam_pass" "$numerator" "$denominator" "$phase"
                done
            done
        done
    done
done
