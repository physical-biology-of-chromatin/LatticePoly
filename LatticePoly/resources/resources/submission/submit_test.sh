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
QUEUE="Lake"

# Max. memory per task
MAXMEM="1G"

# Associated scratch directory
SCRATCHDIR=/home

# Script (relative) path
SCRIPTDIR=$(dirname "$0")

CURRENTDIR=`pwd`

# sbatch arguments
QARGS="--ntasks=1 --mem=${MAXMEM} -J $1 -a 1-$2 -t ${WTIME} -p ${QUEUE}"

# sbatch variables
QVARS="SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR},CURRENTDIR=${CURRENTDIR}"

# Check input parameters and submit
if [ "$#" -eq "2" ]; then
	sbatch ${QARGS} --export=${QVARS} ${SCRIPTDIR}/test.sh
else
	echo -e "\033[1;31mUsage is $0 jobName numJob\033[0m"
fi
