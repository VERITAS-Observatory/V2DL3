#!/bin/bash
#
# Expected files in tar ball:
# - 64080.anasum.root
# - EffectiveArea.root
# - 64080.db.fits.gz
# - ED-pointlike-CI.fits.gz
# - ED-pointlike-CI.fits.gz
# - ED-pointlike-db-CI.fits.gz
# - ED-pointlike-all-CI.fits.gz
# - ED-full-enclosure-all-CI.fits.gz
# - eventfilter.yml
# (plus corresponding log files)

set -e

wget --no-verbose https://syncandshare.desy.de/index.php/s/9P3eaCbqZf7SrdK/download/github-CI.tar.gz
tar -xvzf github-CI.tar.gz
