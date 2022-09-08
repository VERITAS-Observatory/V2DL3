#!/bin/zsh
#
#$ -S /bin/zsh
#
#(Name of the process)
#$ -N v2dl3
#
#$ -l h_cpu=0:20:00
#
#$ -l h_vmem=500M
#
#$ -l tmpdir_size=1G
#
#$ -P cta
#$ -js 9
#
# -j
#
#(use scientific linux 7)
#$ -l os=sl7

COMMAND=$1
CONDA_EXE=$2
ENV=$3
ROOTSYS=$4

# set up root
export PATH=${ROOTSYS}/bin:${PATH}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}
cd $ROOTSYS
source ./bin/thisroot.sh
cd -

# set up conda environment
CONDADIRNAME=$(dirname "$CONDA_EXE")
source $CONDADIRNAME/../etc/profile.d/conda.sh
conda activate $ENV

echo $COMMAND
$COMMAND