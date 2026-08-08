#!/usr/bin/bash

run_types=(hmsdis)
beam_passes=(3 4 5)
targets=(al c cu ld2 lh2 dummy dummy_up dummy_down)
phases=(i ii)

script_path="./MAKE_csvs.py"

rm -f AL/* C/* CU/* DUMMY/* DUMMY_DOWN/* DUMMY_UP/* LD2/* LH2/*

for run_type in "${run_types[@]}"; do
    for beam_pass in "${beam_passes[@]}"; do
        for target in "${targets[@]}"; do
            for phase in "${phases[@]}"; do
                echo "Running $run_type $beam_pass pass $target phase $phase"
                "$script_path" "$run_type" "$beam_pass" "$target" "$phase"
                echo "************************************************"
            done
        done
    done
done

