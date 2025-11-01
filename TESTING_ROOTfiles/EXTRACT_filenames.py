#!/usr/bin/env python3

import os
import re

directory = "/work/hallc/c-rsidis/skimfiles/pass0"
output_file = "filenames.txt"

patterns = [
    re.compile(r"(skimmed_hms_coin_replay_production_\d+_-1\.root)"),
    re.compile(r"(skimmed_shms_coin_replay_production_\d+_-1\.root)"),
    re.compile(r"(skimmed_coin_replay_production_\d+_-1\.root)")
]

def extract_filenames(dirpath):
    filenames = set()
    for filename in os.listdir(dirpath):
        for pattern in patterns:
            match = pattern.match(filename)
            if match:
                filenames.add(match.group(1))
                break
    return sorted(filenames)

def write_filenames(filenames, filepath):
    with open(filepath, "w") as f:
        f.write(",".join(filenames))

def main():
    files = extract_filenames(directory)
    write_filenames(files, output_file)
    print(f"Wrote {len(files)} filenames to {output_file}")

if __name__ == "__main__":
    main()
