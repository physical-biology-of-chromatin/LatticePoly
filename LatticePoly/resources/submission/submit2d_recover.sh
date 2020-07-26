##
##  submit2d_psmn32.sh
##  LatticePoly
##
##  Created by mtortora on 27/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

#!/bin/bash

# Max. walltime
WTIME=168:00:00

# Available queues
PSMN_Q32="CLG6242deb384A,CLG6242deb384B,CLG6242deb384C,\
SLG6142deb384A,SLG6142deb384B,SLG6142deb384C,SLG6142deb384D"

# Associated scratch directory
SCRATCHDIR=/home

# Script (relative) path
SCRIPTDIR=$(dirname "$0")

RECOVERFILE=$1

PARAM1=$(head -n 1 ${RECOVERFILE} | awk '{print $2}')
PARAM2=$(head -n 1 ${RECOVERFILE} | awk '{print $1}')

NJOBS=$(tail -n +2 ${RECOVERFILE} | wc -l)

# qsub arguments
QARGS="-t 1-${NJOBS} -q ${PSMN_Q32} -l h_rt=${WTIME}"

# qsub variables
QVARS="PARAM1=${PARAM1},PARAM2=${PARAM2},RECOVERFILE=${RECOVERFILE},SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR}"

# Check input parameters and submit
if [ "$#" -eq "1" ]; then
	qsub ${QARGS} -N ${PARAM1}${PARAM2}R -v ${QVARS} ${SCRIPTDIR}/sge_recover.sh
else
	echo -e "\033[1;31mUsage is $0 recoverFile\033[0m"
fi
