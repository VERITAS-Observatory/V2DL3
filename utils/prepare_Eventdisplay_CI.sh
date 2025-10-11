#!/bin/bash
# Prepare Eventdisplay CI files for integration tests
#
# kNNInterpolator results are OS depending. Generate
# the files below on a linux machine.

PDIR=$(pwd)
RUNID=64080

echo "Preparing test data in ./test-eventdisplay-CI"
mkdir -p test-eventdisplay-CI
cd test-eventdisplay-CI || exit
wget https://syncandshare.desy.de/index.php/s/9P3eaCbqZf7SrdK/download/github-CI.tar.gz
tar -xvzf github-CI.tar.gz
rm -f ED-*.fits.gz ED-*.log *.tar.gz

# point-like tests
v2dl3-eventdisplay \
    -f ./${RUNID}.anasum.root \
        ./effectiveArea.root \
        --logfile ED-pointlike-CI.log \
        ED-pointlike-CI.fits.gz

# full-enclosure tests
v2dl3-eventdisplay --full-enclosure \
    -f ./${RUNID}.anasum.root \
       ./effectiveArea.root \
       --logfile ED-fullenclosure-CI.log \
       ED-fullenclosure-CI.fits.gz

# point-like tests with DB
v2dl3-eventdisplay \
    -f ./${RUNID}.anasum.root \
        ./effectiveArea.root \
        --db_fits_file ./${RUNID}.db.fits.gz \
        --logfile ED-pointlike-db-CI.log \
        ED-pointlike-db-CI.fits.gz

# point-like tests with event filter
v2dl3-eventdisplay \
    -f ./${RUNID}.allevents.anasum.root \
       ./effectiveArea.root \
       --evt_filter eventfilter.yml \
       --logfile ED-pointlike-all-CI.log \
       ED-pointlike-all-CI.fits.gz

# full-enclosure tests with event selection
v2dl3-eventdisplay --full-enclosure \
     -f ./${RUNID}.allevents.anasum.root \
        ./effectiveArea.root \
        --evt_filter eventfilter.yml \
        --logfile ED-full-enclosure-all-CI.log \
        ED-full-enclosure-all-CI.fits.gz

# Only create archive if at least one file exists
if ls *.root *.gz *.log *.yml 1> /dev/null 2>&1; then
    tar -cvzf github-CI.tar.gz *.root *.gz *.log *.yml
else
    echo "Warning: No files found to archive. Archive not created." >&2
fi

cd ${PDIR} || exit
