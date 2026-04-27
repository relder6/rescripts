#!/usr/bin/bash

run_types=(hmsdis)
beam_passes=(4 5)
targets=(al c cu ld2 lh2)

script_path="./GRAB_model.py"

for run_type in "${run_types[@]}"; do
  for beam_pass in "${beam_passes[@]}"; do
    for target in "${targets[@]}"; do
      "$script_path" "$run_type" "$beam_pass" "$target"
    done
  done
done
