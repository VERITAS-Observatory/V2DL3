#!/bin/bash

# This script should be called from $GITHUB_WORKSPACE (actions/checkout@v3) mounted within a Docker
# image which has been built from the recipe in this directory.

set -e
. ~/py38/bin/activate
. /software/ROOT_build/bin/thisroot.sh
pip install .
throwerror

v2dl3-vegas -h &&\
git config --global --add safe.directory /V2DL3 &&\
git checkout -t origin/main && pip install . &&\
v2dl3-vegas -h \