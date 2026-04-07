#!/bin/bash

##
##  submit2d_slurm.sh
##  LatticePoly
##
##  Created by ppuel on 28/10/2024.
##  Copyright © 2024 ENS Lyon. All rights reserved.
##

# Max. walltime
WTIME=8-00:00:00

# Partition
QUEUE="Cascade"

# Max. memory per task
MAXMEM="1G"

TEMPORARYDIR=/scratch/Cascade/${LOGNAME}/tmp

# Associated scratch directory
SCRATCHDIR=/scratch/Cascade/${LOGNAME}/data

# Script (relative) path
SCRIPTDIR=$(dirname "$0")

XNFSDIR=/Xnfs/lbmcdb/Jost_team/${LOGNAME}

# sbatch arguments
QARGS="--ntasks=1 --mem=${MAXMEM} -a 1-48 -t ${WTIME} -p ${QUEUE}"

# sbatch variables
QVARS="OUTPUTDIR=$1,PARAM=$2,MIN_VAL=$3,MAX_VAL=$4,PARAM2=$5,MIN_VAL2=$6,MAX_VAL2=$7,SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR},TEMPORARYDIR=${TEMPORARYDIR},XNFSDIR=${XNFSDIR}"

EVARRAY=("0.0" "3.0" "100.0")
LDENSARRAY=("0.001" "0.002" "0.005" "0.001")
VALJLLARRAY=("1" "2" "4" "6" "8" "6" "2" "100")
VALJLPARRAY=("1" "2" "4" "6" "8" "2" "6" "100")
VALJLPPARRAY=("1" "2" "4" "6" "8" "2" "6" "100")

# Check input parameters and submit
if [ "$#" -eq "7" ]; then

	# QVARS2="PARAM2=$5,MIN_VAL2=$6,MAX_VAL2=$7,STEP=1,ITERATION=$i"
	# sbatch ${QARGS} -J $2${5}11 --export=${QVARS},${QVARS2} ${SCRIPTDIR}/slurm_sweepV4.sh
	# QVARS2="PARAM2=$5,MIN_VAL2=$6,MAX_VAL2=$7,STEP=2,ITERATION=$i"
	# sbatch ${QARGS} -J $2${5}12 --export=${QVARS},${QVARS2} ${SCRIPTDIR}/slurm_sweepV4.sh
	# QVARS2="PARAM2=$5,MIN_VAL2=$6,MAX_VAL2=$7,STEP=3,ITERATION=$i"
	# sbatch ${QARGS} -J $2${5}13 --export=${QVARS},${QVARS2} ${SCRIPTDIR}/slurm_sweepV4.sh
 # -d afterany:$2$5$($i-1)$s
	for i in $(seq 0 2); do #EV
		for j in $(seq 0 3); do #Ldens
			for k in $(seq 0 7); do #Val
				for s in $(seq 3); do #Jlp Jlpp
					BEGINTIME=$((15*(((($i+1)*4+($j+1))*8+($k+1))*3+$s-124)+1))

					QVARS2="STEP=$s,EV=${EVARRAY[$i]},LDENS=${LDENSARRAY[$j]},VALJLL=${VALJLLARRAY[$k]},VALJLP=${VALJLPARRAY[$k]},VALJLPP=${VALJLPPARRAY[$k]}"
					sbatch ${QARGS} -J $2$5$i$s --begin=now+${BEGINTIME}minutes --export=${QVARS},${QVARS2} ${SCRIPTDIR}/slurm_sweepV5.sh
				done
			done
		done
	done
else
	echo -e "\033[1;31mUsage is $0 OUTPUTDIR paramName1 minVal1 maxVal1 paramName2 minVal2 maxVal2\033[0m"
fi
