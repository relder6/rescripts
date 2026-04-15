#!/usr/bin/bash

module load hcana || { echo "ERROR: failed to load hcana module"; exit 1; }

for window_min in -13200.0; do
    for window_max in -10200.0; do
        echo "window = ($window_min, $window_max)"
        cd /work/hallc/c-rsidis/relder/hallc_replay_rsidis || { echo "ERROR: cannot cd into working dir"; exit 1; }

        sed -i "11s/.*/pdc_tdc_min_win = $window_min, $window_min, $window_min, $window_min, $window_min, $window_min/" PARAM/SHMS/DC/pdc_cuts.param
        sed -i "12s/.*/                  $window_min, $window_min, $window_min, $window_min, $window_min, $window_min/" PARAM/SHMS/DC/pdc_cuts.param
        sed -i "13s/.*/pdc_tdc_max_win = $window_max, $window_max, $window_max, $window_max, $window_max, $window_max/" PARAM/SHMS/DC/pdc_cuts.param
        sed -i "14s/.*/                  $window_max, $window_max, $window_max, $window_max, $window_max, $window_max/" PARAM/SHMS/DC/pdc_cuts.param

        # run the replay
        hcana -l -q 'SCRIPTS/SHMS/PRODUCTION/replay_production_shms_coin.C(25153,50000)'

        # move rootfiles
        cd /volatile/hallc/c-rsidis/relder/ROOTfiles || { echo "cd failed"; exit 1; }
        mv shms_coin_replay_production_25153_50000.root \
            CALIB/shms_coin_replay_production_twmin_${window_min}_twmax${window_max}_25153_50000.root

        echo "*************************"
    done
done
echo "completed; exiting..."
