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



#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/ReplicationAnalysis.py  ${CURRENTDIR}  20000
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_cool.py ${CURRENTDIR}   5010 5 chrIV
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_cool.py ${CURRENTDIR}   5010 5 chrIV
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_all_cool.py ${CURRENTDIR}   5010 2 chrIV
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_all_cool.py ${CURRENTDIR}   5010 5 chrIV
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_all_cool.py ${CURRENTDIR}   5010 3 chrIV
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_all_cool.py ${CURRENTDIR}   5010 4 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_all_cool.py ${CURRENTDIR}   5010 5 chrIV

