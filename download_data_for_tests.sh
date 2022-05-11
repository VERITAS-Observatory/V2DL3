#!/bin/bash

set -e

VERSION="0.3.0"

wget --no-verbose https://syncandshare.desy.de/index.php/s/jYm3WQfp47Yqnyk/download/64080.anasum.root
wget --no-verbose https://syncandshare.desy.de/index.php/s/L3TxQG2bMeA7Ai4/download/effArea-v486-auxv01-CARE_June2020-Cut-NTel2-PointSource-Moderate-TMVA-BDT-GEO-V6_2012_2013a-ATM62-T1234-testCI.root
mv -f effArea-v486-auxv01-CARE_June2020-Cut-NTel2-PointSource-Moderate-TMVA-BDT-GEO-V6_2012_2013a-ATM62-T1234-testCI.root effectiveArea.root
wget --no-verbose https://syncandshare.desy.de/index.php/s/nGnfCMendbKyXzb/download/ED-${VERSION}-pointlike-CI.fits.gz
mv -f ED-${VERSION}-pointlike-CI.fits.gz ED-pointlike-CI.fits.gz
wget --no-verbose https://syncandshare.desy.de/index.php/s/GAn5NL39YFQAkkF/download/ED-${VERSION}-fullenclosure-CI.fits.gz
mv -f ED-${VERSION}-fullenclosure-CI.fits.gz ED-fullenclosure-CI.fits.gz
