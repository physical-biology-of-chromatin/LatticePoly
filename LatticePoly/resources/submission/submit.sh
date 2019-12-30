#!/bin/bash

if [ "$#" -eq "4" ]; then
	qsub -N sweep_$1 -v PARAM=$1,MINVAL=$2,MAXVAL=$3,N=$4 psmn.sh

else
	echo -e "\033[1;31mUsage is $0 param_name min_val max_val num_jobs\033[0m"
	exit

fi
