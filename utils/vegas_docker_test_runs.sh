#!/bin/bash

# This script should be called from $GITHUB_WORKSPACE (actions/checkout@v3) mounted within a Docker
# image which has been built from the recipe in this directory.


# These have been hardcoded into the Docker image after building for now because they cannot be made public.
# The runlists are dynamically made, so any stage5/EA files may be given as inputs.
STAGE5_DIR='/v2dl3-inputs/stage5'
EA_POINTLIKE=1'v2dl3-inputs/eas/ea-pl-1.root'
EA_POINTLIKE_2='v2dl3-inputs/eas/ea-pl-2.root'
EA_FULL_ENCLOSURE_1='v2dl3-inputs/eas/ea-fe-1.root'
EA_FULL_ENCLOSURE_2='v2dl3-inputs/eas/ea-fe-2.root'
EA_EVCLASS_1='v2dl3-inputs/eas/ea-fe-1.root'
EA_EVCLASS_2='v2dl3-inputs/eas/ea-fe-1.root'


# Can override the inputs in as script arguments instead. Useful for downloading inputs to the git runner before mounting to Docker.
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


set -e
. ~/py38/bin/activate

set +e
. /software/ROOT_build/bin/thisroot.sh

set -e
echo "Installing v2dl3-vegas with changes..."
pip install .

set +e
v2dl3-vegas -h


# ---------- TEST RUNS -----------

# Point-like
function pointlike() {
    v2dl3-vegas --point-like  
}

# Point-like 2


# Full-enclosure 1


# Full-enclosure 2


# Single event class


# Multi event class


git config --global --add safe.directory /V2DL3 &&\
git checkout -t origin/main && pip install . &&\
v2dl3-vegas -h \