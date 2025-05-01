#!/bin/bash

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode v08_00_00_80 -q e17:prof

export PATH=/gluster/home/ogregory/nuwro/bin:$PATH
export LD_LIBRARY_PATH=/gluster/home/ogregory/nuwro/bin:$LD_LIBRARY_PATH
cd /gluster/home/ogregory/Project

mkdir -p logs_plots

# extend this list here for all runs
root -b -q dalitzScript.C > "logs_plots/console_dalitz.txt"
root -b -q omegaKScript.C > "logs_plots/console_omegaK.txt"
root -b -q pKScript.C > "logs_plots/console_pK.txt"
root -b -q q2Script.C > "logs_plots/console_q2.txt"
root -b -q thetaScript.C > "logs_plots/console_theta.txt"
root -b -q wModelsScript.C > "logs_plots/console_wModels.txt"
root -b -q wScript.C > "logs_plots/console_w.txt"
