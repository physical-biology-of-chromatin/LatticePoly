##
##  sge.sh
##  LatticePoly
##
##  Created by mtortora on 29/09/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
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

# Executable path
EXEC=bin/lat

# Data directory on local disk
DATDIR=data/output

# Name of output directory (based on job name & task ID)
ID=$(echo ${SGE_TASK_ID} | awk '{printf("%04d\n", $1-1)}')

OUTDIR=${JOB_NAME}${ID}

# Temporary directory on scratch
TMPDIR=${SCRATCHDIR}/${LOGNAME}/LatticeData/${OUTDIR}

# Create output directory if necessary
[ ! -d "${TMPDIR}" ] && mkdir -p ${TMPDIR}

# Substitution strings
DIRSUB="s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${TMPDIR} ;|;"

# Copy input configuration file to output directory, substituting paths and parameter values
sed -e "${DIRSUB}" < data/input.cfg > ${TMPDIR}/input.cfg
#sed -e "${DIRSUB}" < data/output/box.vtp > ${TMPDIR}/box.vtp
#sed -e "${DIRSUB}" < data/output/poly100001.vtp > ${TMPDIR}/poly100001.vtp



# Run
./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out

# analysis


python3 resources/MonomerDist_HiC_window_all.py ${TMPDIR}  2000  2100 3
python3 resources/MonomerDist_HiC_window_all.py ${TMPDIR}  2100  2200 3
python3 resources/MonomerDist_HiC_window_all.py ${TMPDIR}  2200  2300 3
python3 resources/MonomerDist_HiC_window_all.py ${TMPDIR}  2300  2400 3
python3 resources/MonomerDist_HiC_window_all.py ${TMPDIR}  2400  2500 3
python3 resources/MonomerDist_HiC_window_all.py ${TMPDIR}  2500  2600 3
python3 resources/MonomerDist_HiC_window_all.py ${TMPDIR}  2600  2700 3
python3 resources/MonomerDist_HiC_window_all.py ${TMPDIR}  2700  2800 3
python3 resources/MonomerDist_HiC_window_all.py ${TMPDIR}  2800  2900 3
python3 resources/MonomerDist_HiC_window_all.py ${TMPDIR}  2900  3000 3



python3 resources/Poly_Rcmdiff_SCs.py ${TMPDIR}  2000
python3 resources/Mixing_during_replication.py ${TMPDIR}  2000 3














# Move SGE output files to data directory
mv ${SGE_O_WORKDIR}/${JOB_NAME}.e${JOB_ID}.${SGE_TASK_ID} ${TMPDIR}
mv ${SGE_O_WORKDIR}/${JOB_NAME}.o${JOB_ID}.${SGE_TASK_ID} ${TMPDIR}

# Create data directory on local disk
[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Archive output files to home directory
tar --transform "s|^|${OUTDIR}/|" -czvf ${DATDIR}/${OUTDIR}.tar.gz -C ${TMPDIR} .
tar -xzf ${DATDIR}/${OUTDIR}.tar.gz -C data/output/October22/segregation/null/10



# Clean scratch
rm -rf ${TMPDIR}
