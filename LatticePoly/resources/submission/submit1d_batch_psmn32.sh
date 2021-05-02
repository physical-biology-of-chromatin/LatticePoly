##
##  submit1d_batch_psmn32.sh
##  LatticePoly
##
##  Created by mtortora on 15/04/2021.
##  Copyright © 2021 ENS Lyon. All rights reserved.
##

#!/bin/bash

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
QARGS="-t 1-$4 -q ${PSMN_Q32} -l h_rt=${WTIME}"

# qsub variables
QVARS="PARAM=$1,MIN_VAL=$2,MAX_VAL=$3,SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR}"

# Check input parameters and submit
if [ "$#" -eq "5" ]; then
	if [ $( expr $4 % 32 ) -eq "0" ]; then
		for i in $(seq $5); do
			QVARS2="ITER=$i"
			qsub ${QARGS} -N $1$i -v ${QVARS},${QVARS2} ${SCRIPTDIR}/sge_sweep.sh
		done
	else
		echo "numJob must be multiple of 32 (got $4)"
	fi
else
	echo -e "\033[1;31mUsage is $0 paramName minVal maxVal numJob numRep\033[0m"
fi
