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
ROOTDIR=${SCRIPTDIR}/..

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTDIR/lib
PYTHONPATH=$PYTHONPATH:$ROOTDIR/resources
PATH=/home/aabdulla/miniconda3/bin:$PATH

export LD_LIBRARY_PATH
export PYTHONPATH
export PATH

# Set working directory to root
cd ${ROOTDIR}

# Run
python3 /home/aabdulla/3D/ConfineTopo/LatticePoly/LatticePoly/az_resources/MonomerDist_HiC.py ${DATA_DIR} 950 0 10240 
python3 /home/aabdulla/3D/ConfineTopo/LatticePoly/LatticePoly/az_resources/MonomerDist_HiC.py ${DATA_DIR} 500 0 10240 
python3 /home/aabdulla/3D/ConfineTopo/LatticePoly/LatticePoly/az_resources/MonomerDist_HiC.py ${DATA_DIR} 100 0 10240 
python3 /home/aabdulla/3D/ConfineTopo/LatticePoly/LatticePoly/az_resources/MonomerDist_HiC.py ${DATA_DIR} 0 0 10240 