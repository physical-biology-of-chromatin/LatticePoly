##
##  sge_sweep.sh
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
cd ${ROOTDIR}

# Set parameter values by linear interpolation between MIN_VAL and MAX_VAL based on task ID
VAL=$(echo ${MIN_VAL} ${MAX_VAL} ${SGE_TASK_FIRST} ${SGE_TASK_LAST} ${SGE_TASK_ID} | awk '{printf("%.4f\n", $1+($2-$1)/($4-$3)*($5-$3))}')

# Executable path
EXEC=bin/lat

# Data directory on local disk
DATDIR=${PARAM}
[ ! -z "${PARAM2}" ] && DATDIR=${PARAM2}/${VAL2}/${DATDIR}

DATDIR=data${ITER}/${DATDIR}

# Ouput directory on scratch
OUTDIR=${PARAM}_${VAL}
[ ! -z "${PARAM2}" ] && OUTDIR=${PARAM2}_${VAL2}_${OUTDIR}

OUTDIR=${SCRATCHDIR}/${LOGNAME}/LatticeData${ITER}/${OUTDIR}

# Create output directory if necessary
[ ! -d "${OUTDIR}" ] && mkdir -p ${OUTDIR}

# Substitution strings
DIRSUB="s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${OUTDIR} ;|;"
VALSUB="s|\(${PARAM}[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${VAL} ;|;"

[ ! -z "${PARAM2}" ] && VAL2SUB="s|\(${PARAM2}[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${VAL2} ;|;"

# Copy input configuration file to output directory, substituting paths and parameter values
sed -e "${DIRSUB}""${VALSUB}""${VAL2SUB}" < data/input.cfg > ${OUTDIR}/input.cfg

# Run
./${EXEC} ${OUTDIR}/input.cfg > ${OUTDIR}/log.out

# Perform post-processing analyses
python3 resources/DistanceMap.py ${OUTDIR} -1 10 >> ${OUTDIR}/process.out
python3 resources/LiqCluster.py ${OUTDIR} -1 >> ${OUTDIR}/process.out
python3 resources/LiqDensity.py ${OUTDIR} -1 >> ${OUTDIR}/process.out
python3 resources/LiqMSD.py ${OUTDIR} -1 >> ${OUTDIR}/process.out
python3 resources/LiqPolyContact.py ${OUTDIR} -1 >> ${OUTDIR}/process.out
python3 resources/PolyGyration.py ${OUTDIR} -1 >> ${OUTDIR}/process.out
python3 resources/PolyMSD.py ${OUTDIR} -1 >> ${OUTDIR}/process.out

# Move SGE output files to data directory
mv ${SGE_O_WORKDIR}/${JOB_NAME}.e${JOB_ID}.${SGE_TASK_ID} ${OUTDIR}
mv ${SGE_O_WORKDIR}/${JOB_NAME}.o${JOB_ID}.${SGE_TASK_ID} ${OUTDIR}

# Create data directory on local disk
[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Archive output files to home directory
tar --transform "s|^|${VAL}/|" -czvf ${DATDIR}/${VAL}.tar.gz -C ${OUTDIR} .

# Clean scratch
rm -rf ${OUTDIR}
