#!/bin/bash

##
##  submit2d_recover.sh
##  LatticePoly
##
##  Created by mtortora on 27/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

module load Python/3.6.1

source ${HOME}/software/vpython/bin/activate

# Max. walltime
WTIME=168:00:00

# Available queues
PSMN_Q32="Epyc7702deb512,\
CLG6242deb384A,CLG6242deb384B,CLG6242deb384C,\
CLG5218deb192A,CLG5218deb192B,CLG5218deb192C,CLG5218deb192D,\
CLG6226Rdeb192A,CLG6226Rdeb192B,CLG6226Rdeb192C,CLG6226Rdeb192D,\
SLG6142deb384A,SLG6142deb384B,SLG6142deb384C"

# Associated scratch directory
SCRATCHDIR=/scratch/Lake

# Script (relative) path
SCRIPTDIR=$(dirname "$0")

# qsub arguments
QARGS="-q ${PSMN_Q32} -l h_rt=${WTIME}"

# qsub variables
QVARS="SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR}"

# Check input parameters and submit
if [ "$#" -eq "1" ]; then
	RECOVERFILE=$1
	PARAM1=$(head -n 1 ${RECOVERFILE} | awk '{print $2}')
	PARAM2=$(head -n 1 ${RECOVERFILE} | awk '{print $1}')
	NJOBS=$(tail -n +2 ${RECOVERFILE} | wc -l)
	QVARS2="PARAM1=${PARAM1},PARAM2=${PARAM2},RECOVERFILE=${RECOVERFILE}"
	qsub ${QARGS} -t 1-${NJOBS} -N ${PARAM1}${PARAM2}R -v ${QVARS},${QVARS2} ${SCRIPTDIR}/sge_recover.sh
else
	echo -e "\033[1;31mUsage is $0 recoverFile\033[0m"
fi
