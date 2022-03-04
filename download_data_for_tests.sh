#!/bin/bash

set -e

VERSION="0.2.1"

wget https://desycloud.desy.de/index.php/s/fDrSbSYjB4SJ9np/download/64080.anasum.root
wget https://desycloud.desy.de/index.php/s/fEDqTCRgmTbomiG/download/effArea-v486-auxv01-CARE_June2020-Cut-NTel2-PointSource-Moderate-TMVA-BDT-GEO-V6_2012_2013a-ATM62-T1234-testCI.root
mv -f effArea-v486-auxv01-CARE_June2020-Cut-NTel2-PointSource-Moderate-TMVA-BDT-GEO-V6_2012_2013a-ATM62-T1234-testCI.root effectiveArea.root
wget https://desycloud.desy.de/index.php/s/GMRAwjJ2BQnM2LN/download/ED-${VERSION}-pointlike-CI.fits.gz
mv -f ED-${VERSION}-pointlike-CI.fits.gz ED-pointlike-CI.fits.gz
wget https://desycloud.desy.de/index.php/s/rQ3QtEMGsCZbPaw/download/ED-${VERSION}-fullenclosure-CI.fits.gz
mv -f ED-${VERSION}-fullenclosure-CI.fits.gz ED-fullenclosure-CI.fits.gz
