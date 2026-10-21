#!/usr/bin/bash

run_types=(hmsdis)
beam_passes=(3 4 5)
targets=(al c cu dummy dummy_up dummy_down)
targets_cryo=(ld2 lh2)
phases=(i ii)

script_path="./MAKE_csvs.py"

max_jobs=6

rm -rf AL/* C/* CU/* DUMMY/* DUMMY_DOWN/* DUMMY_UP/* LD2/* LH2/*

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for target in "${targets[@]}"; do
            for phase in "${phases[@]}"; do

                if [[ "$beam_pass" == "3" && "$phase" == "i" ]]; then
                    continue
                fi
                
                echo "Running $run_type $beam_pass pass $target phase $phase"
                "$script_path" "$run_type" "$beam_pass" "$target" "$phase" &

                while (( $(jobs -rp | wc -l) >= max_jobs )); do
                    wait -n
                done
            done
        done
    done
done

wait

echo "Solid target jobs done."
echo "Starting cryo target jobs..."

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for target in "${targets_cryo[@]}"; do
            for phase in "${phases[@]}"; do

                if [[ "$beam_pass" == "3" && "$phase" == "i" ]]; then
                    continue
                fi
                
                echo "Running $run_type $beam_pass pass $target phase $phase"
                "$script_path" "$run_type" "$beam_pass" "$target" "$phase" &

                while (( $(jobs -rp | wc -l) >= max_jobs )); do
                    wait -n
                done
            done
        done
    done
done
