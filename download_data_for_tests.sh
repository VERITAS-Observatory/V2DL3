#!/bin/bash

set -e

VERSION="0.3.0"

wget --no-verbose https://syncandshare.desy.de/index.php/s/jYm3WQfp47Yqnyk/download/64080.anasum.root
wget --no-verbose https://syncandshare.desy.de/index.php/s/BTrixGzZAPza8Wb/download/effectiveArea.root
wget --no-verbose https://syncandshare.desy.de/index.php/s/nGnfCMendbKyXzb/download/ED-${VERSION}-pointlike-CI.fits.gz
mv -f ED-${VERSION}-pointlike-CI.fits.gz ED-pointlike-CI.fits.gz
wget --no-verbose https://syncandshare.desy.de/index.php/s/GAn5NL39YFQAkkF/download/ED-${VERSION}-fullenclosure-CI.fits.gz
mv -f ED-${VERSION}-fullenclosure-CI.fits.gz ED-fullenclosure-CI.fits.gz
