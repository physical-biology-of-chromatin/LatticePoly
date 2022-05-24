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
sed -e "${DIRSUB}" < data/box.vtp > ${TMPDIR}/box.vtp
sed -e "${DIRSUB}" < data/poly100001.vtp > ${TMPDIR}/poly100001.vtp



# Run
./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out

# analysis
python3 resources/Yeast_distance_1.py ${TMPDIR}  100101
python3 resources/Yeast_distance_2.py ${TMPDIR}  100101
python3 resources/Yeast_distance_3.py ${TMPDIR}  100101
python3 resources/Forksnumber.py ${TMPDIR}  100101
python3 resources/MonomerDist_HiC.py ${TMPDIR}  100201 0 872   #after G1
python3 resources/MonomerDist_HiC.py ${TMPDIR}  100334 0 872   #1min
python3 resources/MonomerDist_HiC.py${TMPDIR}  100401 0 872   #3min
python3 resources/MonomerDist_HiC.py ${TMPDIR}  100501 0 872   #6min
python3 resources/MonomerDist_HiC.py ${TMPDIR}  100601 0 872  #9min
python3 resources/MonomerDist_HiC.py ${TMPDIR}  100801 0 872   #15min
python3 resources/MonomerDist_HiC.py ${TMPDIR}  101101 0 872   #24min



#python3 resources/Forksnumber.py ${TMPDIR}  100801


# Move SGE output files to data directory
mv ${SGE_O_WORKDIR}/${JOB_NAME}.e${JOB_ID}.${SGE_TASK_ID} ${TMPDIR}
mv ${SGE_O_WORKDIR}/${JOB_NAME}.o${JOB_ID}.${SGE_TASK_ID} ${TMPDIR}

# Create data directory on local disk
[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Archive output files to home directory
tar --transform "s|^|${OUTDIR}/|" -czvf ${DATDIR}/${OUTDIR}.tar.gz -C ${TMPDIR} .

# Clean scratch
rm -rf ${TMPDIR}
