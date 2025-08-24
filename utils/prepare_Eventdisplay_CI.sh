#!/bin/bash
# Prepare Eventdisplay CI files for integration tests
#
# kNNInterpolator results are OS depending. Generate
# the files below on a linux machine.

PDIR=$(pwd)

echo "Preparing test data in ./test-eventdisplay-CI"
mkdir -p test-eventdisplay-CI
cd test-eventdisplay-CI || exit
../download_Eventdisplay_test_data.sh
rm -f ED-*.fits.gz ED-*.log *.tar.gz

# point-like tests
python ../../pyV2DL3/script/v2dl3_for_Eventdisplay.py \
    -f ./64080.anasum.root \
        ./effectiveArea.root \
        --logfile ED-pointlike-CI.log \
        ED-pointlike-CI.fits.gz

# full-enclosure tests
python ../../pyV2DL3/script/v2dl3_for_Eventdisplay.py --full-enclosure \
    -f ./64080.anasum.root \
       ./effectiveArea.root \
       --logfile ED-fullenclosure-CI.log \
       ED-fullenclosure-CI.fits.gz

# point-like tests with DB
python ../../pyV2DL3/script/v2dl3_for_Eventdisplay.py \
    -f ./64080.anasum.root \
        ./effectiveArea.root \
        --db_fits_file ./64080.db.fits.gz \
        --logfile ED-pointlike-db-CI.log \
        ED-pointlike-db-CI.fits.gz

# point-like tests with event filter
python ../../pyV2DL3/script/v2dl3_for_Eventdisplay.py \
    -f ./64080.allevents.anasum.root \
       ./effectiveArea.root \
       --evt_filter eventfilter.yml \
       --logfile ED-pointlike-all-CI.log \
       ED-pointlike-all-CI.fits.gz

# full-enclosure tests with event selection
python ../../pyV2DL3/script/v2dl3_for_Eventdisplay.py --full-enclosure \
     -f ./64080.allevents.anasum.root \
        ./effectiveArea.root \
        --evt_filter eventfilter.yml \
        --logfile ED-full-enclosure-all-CI.log \
        ED-full-enclosure-all-CI.fits.gz

tar -cvzf github-CI.tar.gz *.root *.gz *.log *.yml

cd ${PDIR} || exit
