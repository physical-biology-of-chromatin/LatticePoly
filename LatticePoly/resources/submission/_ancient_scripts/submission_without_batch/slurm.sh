#!/bin/bash

##
##  slurm.sh
##  LatticePoly
##
##  Created by mtortora on 20/06/2022
##  Copyright © 2022 ENS Lyon. All rights reserved.
##

# Set working directory to root
cd ~/Documents/LatticePoly/LatticePoly

# Executable path
EXEC=bin/lat

# Data directory on local disk
DATDIR=data/output/testVTK/

# Name of output directory (based on job name & task ID)
ID=1

OUTDIR=test${ID}

# Temporary directory on scratch
TMPDIR=${SCRATCHDIR}/${LOGNAME}/LatticeData/${OUTDIR}

# Create output directory if necessary
[ ! -d "${TMPDIR}" ] && mkdir -p ${TMPDIR}

# Substitution strings
DIRSUB="s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${TMPDIR} ;|;"

# Copy input configuration file to output directory, substituting paths and parameter values
sed -e "${DIRSUB}" < data/input.cfg > ${TMPDIR}/input.cfg

# Run
./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out

# Move SLURM output files to data directory
# mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${TMPDIR}
# mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${TMPDIR}

# Create data directory on local disk
[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Archive output files to home directory
tar --transform "s|^|${OUTDIR}/|" -czf ${DATDIR}/${OUTDIR}.tar.gz -C ${TMPDIR} .

# Clean scratch
rm -rf ${TMPDIR}
