#!/bin/bash

apptainer run \
          -B /cvmfs \
          -B /run/user \
          -B /etc/hosts \
          -B /etc/localtime \
    -B /gluster/data \
          --env UPS_OVERRIDE='-H Linux64bit+3.10-2.17' \
          /cvmfs/uboone.opensciencegrid.org/containers/uboone-devel-sl7 \
          /gluster/home/ogregory/scripts/run_nuwro.sh
