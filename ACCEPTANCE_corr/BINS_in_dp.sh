#!/usr/bin/bash

for dp_min in -8.0 -7.0 -6.0 -5.0 -4.0 -3.0 -2.0 -1.0 0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0; do
    cd /w/hallc-scshelf2102/c-rsidis/relder/rescripts/INIT || { echo "ERROR: cannot cd into working dir"; exit 1; }

    dp_max=$(echo "$dp_min + 1.0" | bc)
    
    echo "dp range = $dp_min, $dp_max"

    # update delta range in config file

    sed -i '10s/.*/            "USING_DELTA_CORR": False,}/' config.py
    sed -i "16s/.*/    return {\"H_gtr_dp_min_cut\": $dp_min,/" config.py
    sed -i "17s/.*/            \"H_gtr_dp_max_cut\": $dp_max,/" config.py

    cd /w/hallc-scshelf2102/c-rsidis/relder/rescripts/XSEC/MAKE_csvs || { echo "ERROR: cannot cd into working dir"; exit 1; }

    rm -rf AL C CU DUMMY DUMMY_DOWN DUMMY_UP LD2 LH2
    mkdir AL C CU DUMMY DUMMY_DOWN DUMMY_UP LD2 LH2
    ./test_2.sh

    cd /w/hallc-scshelf2102/c-rsidis/relder/rescripts/XSEC/DATA_to_MC || { echo "ERROR: cannot cd into working dir"; exit 1; }

    rm -rf AL C CU DUMMY DUMMY_DOWN DUMMY_UP LD2 LH2 PDFs
    mkdir AL C CU DUMMY DUMMY_DOWN DUMMY_UP LD2 LH2 PDFs
    ./test_2.sh

    # move the data_to_mc files to acceptance_corr directory
    for folder in AL C CU DUMMY DUMMY_DOWN DUMMY_UP LD2 LH2 PDFs; do
        rm -rf "/w/hallc-scshelf2102/c-rsidis/relder/rescripts/ACCEPTANCE_corr/DATA_to_MC/${folder}_${dp_min}_to_${dp_max}"
        mv "$folder" "/w/hallc-scshelf2102/c-rsidis/relder/rescripts/ACCEPTANCE_corr/DATA_to_MC/${folder}_${dp_min}_to_${dp_max}"
        mkdir "$folder"
    done

    # move the make_csvs files to acceptance_corr directory
    cd /w/hallc-scshelf2102/c-rsidis/relder/rescripts/XSEC/MAKE_csvs || { echo "ERROR: cannot cd into working dir"; exit 1; }

    for folder in AL C CU DUMMY DUMMY_DOWN DUMMY_UP LD2 LH2; do
        rm -rf "/w/hallc-scshelf2102/c-rsidis/relder/rescripts/ACCEPTANCE_corr/CSVs/${folder}_${dp_min}_to_${dp_max}"
        mv "$folder" "/w/hallc-scshelf2102/c-rsidis/relder/rescripts/ACCEPTANCE_corr/CSVs/${folder}_${dp_min}_to_${dp_max}"
        mkdir "$folder"
    done

done

cd /w/hallc-scshelf2102/c-rsidis/relder/rescripts/INIT || { echo "ERROR: cannot cd into working dir"; exit 1; }

sed -i '10s/.*/            "USING_DELTA_CORR": True,}/' config.py
sed -i "16s/.*/    return {\"H_gtr_dp_min_cut\": -8.0,/" config.py
sed -i "17s/.*/            \"H_gtr_dp_max_cut\": 8.0,/" config.py
