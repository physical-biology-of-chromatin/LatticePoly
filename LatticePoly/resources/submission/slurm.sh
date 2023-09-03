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

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTDIR/lib
export LD_LIBRARY_PATH

# Set working directory to root
cd ${ROOTDIR}

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTDIR/lib
PYTHONPATH=$PYTHONPATH:$ROOTDIR/resources

export LD_LIBRARY_PATH
export PYTHONPATH

# Executable path
EXEC=bin/lat

# Temporary directory on tmp
TMPDIR=${TEMPORARYDIR}/${OUTPUTDIR}_${SLURM_ARRAY_TASK_ID}

# Create temporary directory if necessary
[ ! -d "${TMPDIR}" ] && mkdir -p ${TMPDIR}

# Scratch directory on scratch
SCRDIR=${SCRATCHDIR}/${OUTPUTDIR}/${SLURM_ARRAY_TASK_ID}

# Create scratch directory if necessary
[ ! -d "${SCRDIR}" ] && mkdir -p ${SCRDIR}

# Data directory on Xnfs
DATDIR=${XNFSDIR}/${OUTPUTDIR}

# Create data directory if necessary
[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Substitution strings
DIRSUB="s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${TMPDIR} ;|;"

# Copy input configuration file to output directory, substituting paths and parameter values
sed -e "${DIRSUB}" < data/input.cfg > ${TMPDIR}/input.cfg

# Run
./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out

# Perform post-processing analyses
python3 resources/DistanceMap.py ${TMPDIR} -1 10 >> ${TMPDIR}/process.out
python3 resources/LiqDensity.py ${TMPDIR} -1 >> ${TMPDIR}/process.out
python3 resources/LiqCluster.py ${TMPDIR} -1 >> ${TMPDIR}/process.out
python3 resources/LiqMSD.py ${TMPDIR} -1 >> ${TMPDIR}/process.out
python3 resources/LiqPolyCoincidence.py ${TMPDIR} -1 >> ${TMPDIR}/process.out
python3 resources/PolyGyration.py ${TMPDIR} -1 >> ${TMPDIR}/process.out
python3 resources/PolyMSD.py ${TMPDIR} -1 >> ${TMPDIR}/process.out
python3 resources/LiqCluster_Coating.py ${TMPDIR} -1 >> ${TMPDIR}/process.out

# Move SLURM output files to data directory
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${TMPDIR}
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${TMPDIR}

# Copy post-processing files to scratch for averaging
cp ${TMPDIR}/liqDropNum.res ${SCRDIR}
cp ${TMPDIR}/liqDropRad.res ${SCRDIR}
cp ${TMPDIR}/liqRadii.res ${SCRDIR}
cp ${TMPDIR}/liqFraction.res ${SCRDIR}
cp ${TMPDIR}/liqMSD.res ${SCRDIR}
cp ${TMPDIR}/liq_Het_Coincidence.res ${SCRDIR}
cp ${TMPDIR}/poly_Het_Coincidence.res ${SCRDIR}
cp ${TMPDIR}/liq_PRE_Coincidence.res ${SCRDIR}
cp ${TMPDIR}/poly_PRE_Coincidence.res ${SCRDIR}
cp ${TMPDIR}/polyAniso.res ${SCRDIR}
cp ${TMPDIR}/polyGyration.res ${SCRDIR}
cp ${TMPDIR}/polyHetMSD.res ${SCRDIR}
cp ${TMPDIR}/polyHomMSD.res ${SCRDIR}
cp ${TMPDIR}/polyPREMSD.res ${SCRDIR}
cp ${TMPDIR}/liqMean.res ${SCRDIR}
cp ${TMPDIR}/liqSTD.res ${SCRDIR}
cp ${TMPDIR}/polyFirstZone.res ${SCRDIR}
cp ${TMPDIR}/polySecondZone.res ${SCRDIR}
cp ${TMPDIR}/polyThirdZone.res ${SCRDIR}

# Archive output files to home directory
tar -czf ${DATDIR}/${SLURM_ARRAY_TASK_ID}.tar.gz ${TMPDIR}/

# Clean scratch
rm -rf ${TMPDIR}
