#!/bin/bash

##
##  submit_slurm.sh
##  LatticePoly
##
##  Created by mtortora on 20/06/2022.
##  Copyright © 2022 ENS Lyon. All rights reserved.
##

# Max. walltime
WTIME=8-00:00:00

# Partition
QUEUE="Cascade"

# Max. memory per task
MAXMEM="1G"

# Associated scratch directory
SCRATCHDIR=/home

# Script (relative) path
SCRIPTDIR=$(dirname "$0")

# sbatch arguments
QARGS="--ntasks=1 --mem=${MAXMEM} -a 1-1 -t ${WTIME} -p ${QUEUE}"

# sbatch variables
#QVARS

# Check input parameters and submit
if [ "$#" -eq "2" ]; then

	for dir in "$2"/*; do
           QVARS="DATA_DIR=${dir},SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR}"
		   sbatch ${QARGS} -J $1 --export=${QVARS} ${SCRIPTDIR}/slurm_analysis.sh
    done

else
	echo -e "\033[1;31mUsage is $0 JobName DataFolder\033[0m"
fi
