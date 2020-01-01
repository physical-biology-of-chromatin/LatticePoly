##
##  sge.sh
##  LatticePoly
##
##  Created by mtortora on 27/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

#$ -S /bin/bash

### Load the user environment for SGE
#$ -cwd

### Export environment variables to all runtime nodes
#$ -V

# Relative path to code root directory
ROOTDIR=${SCRIPTDIR}/../..

# Set working directory to root
cd ${SGE_O_WORKDIR}
cd ${ROOTDIR}

# Set parameter values by linear interpolation between MIN_VAL and MAX_VAL based on task ID
VAL=$(echo ${MIN_VAL} ${MAX_VAL} ${SGE_TASK_FIRST} ${SGE_TASK_LAST} ${SGE_TASK_ID} | awk '{printf("%.3f\n", $1+($2-$1)/($4-$3)*($5-$3))}')

# Executable path
EXEC=bin/lat

# Ouput directory
OUTDIR=${SCRATCHDIR}/${LOGNAME}/LatticeData/${PARAM}_${VAL}

# Create output directory if necessary
if [ ! -d "${OUTDIR}" ]; then
	mkdir -p ${OUTDIR}
fi

# Copy input configuration file to output directory, substituting paths and parameter values
sed -e \
"s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${OUTDIR} ;|;"\
"s|\(${PARAM}[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${VAL} ;|"\
< data/input.cfg > ${OUTDIR}/.input.cfg

# Run
./${EXEC} ${OUTDIR}/.input.cfg > ${OUTDIR}/log.out

# Move SGE output files to data directory
mv ${SGE_O_WORKDIR}/${JOB_NAME}.e${JOB_ID}.${SGE_TASK_ID} ${OUTDIR}
mv ${SGE_O_WORKDIR}/${JOB_NAME}.o${JOB_ID}.${SGE_TASK_ID} ${OUTDIR}

# Create data folder in home directory
if [ ! -d "data/${PARAM}" ]; then
	mkdir -p data/${PARAM}
fi

# Archive output files to home directory
tar --transform "s|^|${VAL}/|" -czvf data/${PARAM}/${VAL}.tar.gz -C ${OUTDIR} .

# Clean scratch
rm -r ${OUTDIR}
