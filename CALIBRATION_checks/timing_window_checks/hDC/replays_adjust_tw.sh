#!/usr/bin/bash

module load hcana || { echo "ERROR: failed to load hcana module"; exit 1; }

for window_min in -13500.0; do
    for window_max in -11400; do
    # for window_max in -11400.0 -11800.0 -11900.0 -12000.0 -12100.0; do
        echo "window = ($window_min, $window_max)"
        cd /work/hallc/c-rsidis/relder/hallc_replay_rsidis || { echo "ERROR: cannot cd into working dir"; exit 1; }

        sed -i "11s/.*/hdc_tdc_min_win = $window_min, $window_min, $window_min, $window_min, $window_min, $window_min/" PARAM/HMS/DC/hdc_cuts.param
        sed -i "12s/.*/                  $window_min, $window_min, $window_min, $window_min, $window_min, $window_min/" PARAM/HMS/DC/hdc_cuts.param
        sed -i "13s/.*/hdc_tdc_max_win = $window_max, $window_max, $window_max, $window_max, $window_max, $window_max/" PARAM/HMS/DC/hdc_cuts.param
        sed -i "14s/.*/                  $window_max, $window_max, $window_max, $window_max, $window_max, $window_max/" PARAM/HMS/DC/hdc_cuts.param

        # run the replay
        hcana -l -q 'SCRIPTS/HMS/PRODUCTION/replay_production_hms_coin.C(24636,100000)'

        # move rootfiles
        cd /volatile/hallc/c-rsidis/relder/ROOTfiles || { echo "cd failed"; exit 1; }
        mv hms_coin_replay_production_24636_100000.root \
            CALIB/hms_coin_replay_production_twmin_${window_min}_twmax${window_max}_24636_100000.root

        echo "*************************"
    done
done
echo "completed; exiting..."
