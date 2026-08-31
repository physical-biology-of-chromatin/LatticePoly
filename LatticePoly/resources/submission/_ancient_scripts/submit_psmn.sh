##
##  submit_psmn.sh
##  LatticePoly
##
##  Created by mtortora on 29/09/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
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
QARGS="-N $1 -t 1-$2 -q ${PSMN_Q32} -l h_rt=${WTIME}"

# qsub variables
QVARS="SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR}"

# Check input parameters and submit
if [ "$#" -eq "2" ]; then
	qsub ${QARGS} -v ${QVARS} ${SCRIPTDIR}/sge.sh
else
	echo -e "\033[1;31mUsage is $0 jobName numJob\033[0m"
fi
