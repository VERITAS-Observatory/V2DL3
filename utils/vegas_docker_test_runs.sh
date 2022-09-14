#!/bin/bash

# This script should be called from $GITHUB_WORKSPACE (actions/checkout@v3) mounted within a Docker
# image which has been built from the recipe in this directory.

# Provide argument "OUTDIR" for output location

# These have been hardcoded into the Docker image after building for now because they cannot be made public.
# The runlists are dynamically made, so any stage5/EA files may be given as inputs instead.
STAGE5_DIR='/v2dl3-inputs/stage5'
EA_POINTLIKE_1='/v2dl3-inputs/eas/ea_tSq_0p01.root'
EA_POINTLIKE_2='/v2dl3-inputs/eas/ea_tSq_0p005.root'
EA_FULL_ENCLOSURE_1='/v2dl3-inputs/eas/ea_tSq129600_MSW0p8to1p1.root'
EA_FULL_ENCLOSURE_2='/v2dl3-inputs/eas/ea_tSq129600_MSW1p1to1p3.root'
EA_EVCLASS_1='/v2dl3-inputs/eas/ea_tSq129600_MSW0p8to1p1.root'
EA_EVCLASS_2='/v2dl3-inputs/eas/ea_tSq129600_MSW1p1to1p3.root'


# Can override the inputs above as named arguments instead. Useful for downloading inputs to the git runner before mounting to Docker.
for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)

   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"

   export "$KEY"="$VALUE"
done

if [ -z "$STAGE5_DIR" ]
then
    echo "must provide STAGE5_DIR"
    exit 0
fi

export LC_ALL=C.UTF-8
export LANG=C.UTF-8
. /software/ROOT_build/bin/thisroot.sh

set -e
echo "Installing v2dl3-vegas with changes..."
pip install . 
set +e

# ---------- TEST RUNS -----------
function run_tests() 
{
    # First func arg is base output dir
    OUTDIR=$1
    EXTRA_FLAGS="-m --save_msw_msl"

    echo "-------------------------------"
    echo "Point-like 1 - Min flags"
    echo "-------------------------------"
    python3 utils/vegas_runlister.py runlist.txt -rd $STAGE5_DIR -e $EA_POINTLIKE_1 --no_prompt
    v2dl3-vegas --point-like -l runlist.txt $OUTDIR/point-like-1

    echo "-------------------------------"
    echo "Point-like 2 - Extra flags"
    echo "-------------------------------"
    python3 utils/vegas_runlister.py runlist.txt -rd $STAGE5_DIR -e $EA_POINTLIKE_2 --no_prompt
    v2dl3-vegas $EXTRA_FLAGS --point-like -l runlist.txt $OUTDIR/point-like-2

    # Coming in King PSF update....

    # echo "-------------------------------"
    # echo "Full-enclosure 1 - Min flags"
    # echo "-------------------------------"
    #python utils/vegas_runlister.py runlist.txt -rd $STAGE5_DIR -e $EA_FULL_ENCLOSURE_1 --no_prompt
    #v2dl3-vegas --full-enclosure -l runlist.txt $OUTDIR/full-enclosure-1

    # echo "-------------------------------"
    # echo "Full-enclosure 2 - Extra flags"
    # echo "-------------------------------"
    #python utils/vegas_runlister.py runlist.txt -rd $STAGE5_DIR -e $EA_FULL_ENCLOSURE_2 --no_prompt
    #v2dl3-vegas $EXTRA_FLAGS --full-enclosure -l runlist.txt $OUTDIR/full-enclosure-2

    echo "-------------------------------"    
    echo "Single event class"
    echo "-------------------------------"

    python3 utils/vegas_runlister.py runlist.txt -rd $STAGE5_DIR -e $EA_EVCLASS_1 --no_prompt
    v2dl3-vegas -ec -l runlist.txt $OUTDIR/single-evclass

    echo "-------------------------------"
    echo "Multi event class - Extra flags"
    echo "-------------------------------"
    python3 utils/vegas_runlister.py runlist.txt -rd $STAGE5_DIR -e $EA_EVCLASS_1 -e $EA_EVCLASS_2 --no_prompt
    v2dl3-vegas $EXTRA_FLAGS -ec -l runlist.txt $OUTDIR/multi-evclass
}

run_tests $OUTDIR

# echo "Installing v2dl3-vegas main branch..."
# git config --global --add safe.directory /V2DL3 &&\
# git checkout -t origin/main && pip install . &&\
# run_tests control