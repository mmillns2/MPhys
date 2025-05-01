#!/bin/bash

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode v08_00_00_80 -q e17:prof

export PATH=/gluster/home/ogregory/nuwro/bin:$PATH
export LD_LIBRARY_PATH=/gluster/home/ogregory/nuwro/bin:$LD_LIBRARY_PATH
cd /gluster/home/ogregory/nuwro

mkdir -p logs_experiments

# extend this list here for all runs
nuwro -i params.txt -o condor_test.root > "logs_experiments/console_params.txt"
nuwro -i params2.txt -o condor_test2.root > "logs_experiments/console_params2.txt"
