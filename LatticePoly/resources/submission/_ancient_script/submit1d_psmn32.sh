##
##  submit1d_psmn32.sh
##  LatticePoly
##
##  Created by mtortora on 27/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

#!/bin/bash

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
QARGS="-t 1-$4 -q ${PSMN_Q32} -l h_rt=${WTIME}"

# qsub variables
QVARS="PARAM=$1,MIN_VAL=$2,MAX_VAL=$3,SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR}"

# Check input parameters and submit
if [ "$#" -eq "4" ]; then
	if [ $( expr $4 % 32 ) -eq "0" ]; then
		qsub ${QARGS} -N $1 -v ${QVARS} ${SCRIPTDIR}/sge_sweep.sh
	else
		echo "numJob must be multiple of 32 (got $4)"
	fi
else
	echo -e "\033[1;31mUsage is $0 paramName minVal maxVal numJob\033[0m"
fi
