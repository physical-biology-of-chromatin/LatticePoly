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
PSMN_Q32="Epyc7702deb512,c6420deb192,c6420deb192,xlr178deb192,xlr178deb192,xlr170deb384,xlr170deb384,c6420deb384,c6420deb384"
# Associated scratch directory
SCRATCHDIR=/home

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
