#!/bin/bash


##
##  slurm.sh
##  LatticePoly
##
##  Created by mtortora on 20/06/2022
##  Copyright © 2022 ENS Lyon. All rights reserved.
##

#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# Relative path to code root directory
ROOTDIR=${SCRIPTDIR}/../..

# Set working directory to root



python3 /home/ddasaro//LatticePoly/LatticePoly/resources/ReplicationAnalysis.py  ${CURRENTDIR}  20000
