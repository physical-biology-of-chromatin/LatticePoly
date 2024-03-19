#!/bin/bash

##
##  slurm.sh
##  LatticePoly
##
##  Created by mtortora on 20/06/2022
##  Copyright © 2022 ENS Lyon. All rights reserved.
##

#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# Relative path to code root directory
ROOTDIR=${SCRIPTDIR}/../..

# Set working directory to root
cd ${ROOTDIR}

# Executable path
EXEC=bin/lat

# Data directory on local disk
DATDIR=data/output

# Name of output directory (based on job name & task ID)
ID=$(echo ${SLURM_ARRAY_TASK_ID} | awk '{printf("%04d\n", $1-1)}')

OUTDIR=${SLURM_JOB_NAME}${ID}

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

#Analysis
python3 resources/MonomerDist_HiC_M_all.py ${TMPDIR}  10000 1
python3 resources/MonomerDist_HiC_M_cis.py ${TMPDIR}  10000 1
python3 resources/MonomerDist_HiC_M_trans.py ${TMPDIR}  10000 1
python3 resources/MonomerDist_HiC_M_all.py ${TMPDIR}  10000 2
python3 resources/MonomerDist_HiC_M_cis.py ${TMPDIR}  10000 2
python3 resources/MonomerDist_HiC_M_trans.py ${TMPDIR}  10000 2
python3 resources/MonomerDist_HiC_M_all.py ${TMPDIR}  10000 3
python3 resources/MonomerDist_HiC_M_cis.py ${TMPDIR}  10000 3
python3 resources/MonomerDist_HiC_M_trans.py ${TMPDIR}  10000 3
python3 resources/MonomerDist_HiC_M_all.py ${TMPDIR}  10000 4
python3 resources/MonomerDist_HiC_M_cis.py ${TMPDIR}  10000 4
python3 resources/MonomerDist_HiC_M_trans.py ${TMPDIR}  10000 4
python3 resources/MonomerDist_HiC_M_all.py ${TMPDIR}  10000 5
python3 resources/MonomerDist_HiC_M_cis.py ${TMPDIR}  10000 5
python3 resources/MonomerDist_HiC_M_trans.py ${TMPDIR}  10000 5

# Move SLURM output files to data directory
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${TMPDIR}
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${TMPDIR}

# Create data directory on local disk
[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Archive output files to home directory
tar --transform "s|^|${OUTDIR}/|" -czf ${DATDIR}/${OUTDIR}.tar.gz -C ${TMPDIR} .
tar -xzf ${DATDIR}/${OUTDIR}.tar.gz -C data/output/December22/cohesion/

# Clean scratch
rm -rf ${TMPDIR}
