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
SCRATCHDIR=/scratch/Cascade/${LOGNAME}/data

# Script (relative) path
SCRIPTDIR=$(dirname "$0")

# sbatch arguments
QARGS="--ntasks=1 --mem=${MAXMEM} -a 1-48 -t ${WTIME} -p ${QUEUE}"

# sbatch variables
QVARS="OUTPUTDIR=$1,PARAM=$2,MIN_VAL=$3,MAX_VAL=$4,SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR}"

# Check input parameters and submit
if [ "$#" -eq "8" ]; then
	for i in $(seq $8); do
		for s in $(seq 3); do
			QVARS2="PARAM2=$5,MIN_VAL2=$6,MAX_VAL2=$7,STEP=$s,ITERATION=$i"
			sbatch ${QARGS} -J $2$5$i$s --export=${QVARS},${QVARS2} ${SCRIPTDIR}/slurm_sweepV3.sh
		done
	done
else
	echo -e "\033[1;31mUsage is $0 OUTPUTDIR paramName1 minVal1 maxVal1 paramName2 minVal2 maxVal2 numRep\033[0m"
fi
