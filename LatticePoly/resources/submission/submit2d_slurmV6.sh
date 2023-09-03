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

GENE="BX"

TEMPORARYDIR=/scratch/Cascade/${LOGNAME}/tmp

# Associated scratch directory
SCRATCHDIR=/scratch/Cascade/${LOGNAME}/data

# Script (relative) path
SCRIPTDIR=$(dirname "$0")

XNFSDIR=/Xnfs/lbmcdb/Jost_team/${LOGNAME}

# sbatch arguments
QARGS="--ntasks=1 --mem=${MAXMEM} -a 1-48 -t ${WTIME} -p ${QUEUE}"

# sbatch variables
QVARS="OUTPUTDIR=$1,SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR},TEMPORARYDIR=${TEMPORARYDIR},XNFSDIR=${XNFSDIR},GENE=${GENE}" #,EVARRAY=${EVARRAY},JLPARRAY=${JLPARRAY}"




# Check input parameters and submit
if [ "$#" -eq "1" ]; then

	for i in $(seq 15); do #FIELD
        for j in $(seq 3); do #EV
            BEGINTIME=$((10*($i-1)+1))
                    
            
            QVARS2="NUMFIELD=$i,NUM_STEP=$j"
            
            sbatch ${QARGS} -J job$i --begin=now+${BEGINTIME}minutes --export=${QVARS},${QVARS2} ${SCRIPTDIR}/slurm_sweepV6.sh

        done
	done

else
	echo -e "\033[1;31mUsage is $0 OUTPUTDIR\033[0m"
fi
