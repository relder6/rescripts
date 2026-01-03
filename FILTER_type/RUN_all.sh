#!/usr/bin/bash

run_types=(hmsdis)
beam_passes=(1 4 5)
targets=(al c cu ld2 lh2 dummy)

script_path="./FILTER_type.py"

for run_type in "${run_types[@]}"; do
  for beam_pass in "${beam_passes[@]}"; do
    for target in "${targets[@]}"; do
      "$script_path" "$run_type" "$beam_pass" "$target"
    done
  done
done

