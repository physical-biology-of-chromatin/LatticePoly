#!/bin/bash

CUR_DIR=`dirname "$0"`

if [ "$#" -eq "4" ]; then
	qsub -N $1 -t 1-$4 -v PARAM=$1,MIN_VAL=$2,MAX_VAL=$3 ${CUR_DIR}/psmn.sh

else
	echo -e "\033[1;31mUsage is $0 param_name min_val max_val num_jobs\033[0m"
	exit

fi
