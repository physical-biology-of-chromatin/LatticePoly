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

# Associated temporary directory
TEMPORARYDIR=/home/${LOGNAME}/tmp

# Associated scratch directory
SCRATCHDIR=/scratch/Cascade/${LOGNAME}/data

# Script (relative) path
SCRIPTDIR=$(dirname "$0")

# Associated data directory
XNFSDIR=/Xnfs/lbmcdb/Jost_team/${LOGNAME}

# sbatch arguments
QARGS="--ntasks=1 --mem=${MAXMEM} -J $1 -a 1-$2 -t ${WTIME} -p ${QUEUE} --mail-type=END --mail-user=paul-swann.puel@ens-lyon.fr"

# sbatch variables
QVARS="OUTPUTDIR=$1,SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR},TEMPORARYDIR=${TEMPORARYDIR},XNFSDIR=${XNFSDIR}"


# Check input parameters and submit
if [ "$#" -eq "2" ]; then
	sbatch ${QARGS} --export=${QVARS} ${SCRIPTDIR}/slurm.sh
else
	echo -e "\033[1;31mUsage is $0 OUTPUTDIR numRep\033[0m"
fi
