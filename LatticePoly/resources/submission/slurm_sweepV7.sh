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


LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTDIR/lib
PYTHONPATH=$PYTHONPATH:$ROOTDIR/resources


export LD_LIBRARY_PATH
export PYTHONPATH


# Executable path
EXEC=bin/lat

#sbatch variables
# QVARS="OUTPUTDIR=$1,SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR},TEMPORARYDIR=${TEMPORARYDIR},XNFSDIR=${XNFSDIR},GENE=${GENE},EVARRAY=${EVARRAY},JLPARRAY=${JLPARRAY}"
# QVARS2="NUMFIELD=$i"

JLPP_INDEX=${SLURM_ARRAY_TASK_ID}

JLPARRAY=("0.0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" "1.6" "1.7" "1.8" "1.9" "2.0" "2.1" "2.2" "2.3" "2.4" "2.5" "2.6" "2.7" "2.8" "2.9" "3.0" "3.1" "3.2" "3.3" "3.4" "3.5" "3.6" "3.7" "3.8" "3.9" "4.0" "4.1" "4.2" "4.3" "4.4" "4.5" "4.6" "4.7")
JLPPARRAY=("0.0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1.0" "1.1" "1.2" "1.3" "1.4" "1.5" "1.6" "1.7" "1.8" "1.9" "2.0" "2.1" "2.2" "2.3" "2.4" "2.5" "2.6" "2.7" "2.8" "2.9" "3.0" "3.1" "3.2" "3.3" "3.4" "3.5" "3.6" "3.7" "3.8" "3.9" "4.0" "4.1" "4.2" "4.3" "4.4" "4.5" "4.6" "4.7")

JLP=${JLPARRAY[${JLP_INDEX}]}
# echo ${JLP} >> ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out

JLPP=${JLPPARRAY[${JLPP_INDEX}]}
# echo ${EV} >> ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out

# Temporary directory on tmp
TMPDIR=JLP_${JLP}_JLPP_${JLPP}_STEP_${STEP}

TMPDIR=${TEMPORARYDIR}/${OUTPUTDIR}_${TMPDIR}

# Create temporary directory if necessary
[ ! -d "${TMPDIR}" ] && mkdir -p ${TMPDIR}

# Scratch directory on scratch
SCRDIR=JLP/${JLP}/JLPP/${JLPP}/STEP/${STEP}

SCRDIR=${SCRATCHDIR}/${OUTPUTDIR}/${SCRDIR}

# Create scratch directory if necessary
[ ! -d "${SCRDIR}" ] && mkdir -p ${SCRDIR}

# Data directory on Xnfs
DATDIR=JLP/${JLP}/JLPINDEX/${JLP_INDEX}/JLPP/${JLPP}/JLPPINDEX/${JLPP_INDEX}/STEP/


DATDIR=${XNFSDIR}/${OUTPUTDIR}/${DATDIR}

# Create data directory if necessary
[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}



# Substitution strings
DIRSUB="s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${TMPDIR} ;|;"
JLPSUB="s|\(Jlp[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JLP} ;|;"

JLPPSUB="s|\(Jlpp[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JLPP} ;|;"



# Copy input configuration file to output directory, substituting paths and parameter values
sed -e "${DIRSUB}""${JLPSUB}""${JLPPSUB}" < data/input.cfg > ${TMPDIR}/input.cfg


# Run
./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out

# Perform post-processing analyses
python3 resources/LiqCluster.py ${TMPDIR} -1 >> ${TMPDIR}/process.out
python3 resources/PolyMSD.py ${TMPDIR} -1 >> ${TMPDIR}/process.out



# Move SLURM output files to data directory
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${TMPDIR}
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${TMPDIR}

# Copy post-processing files to scratch for averaging
# [ -d "${TMPDIR}/liqDropNum.res" ] && cp ${TMPDIR}/liqDropNum.res ${SCRDIR}
# [ -d "${TMPDIR}/liqDropRad.res" ] && cp ${TMPDIR}/liqDropRad.res ${SCRDIR}
# [ -d "${TMPDIR}/liqRadii.res" ] && cp ${TMPDIR}/liqRadii.res ${SCRDIR}
# [ -d "${TMPDIR}/liqFraction.res" ] && cp ${TMPDIR}/liqFraction.res ${SCRDIR}
# [ -d "${TMPDIR}/liqMSD.res" ] && cp ${TMPDIR}/liqMSD.res ${SCRDIR}
# [ -d "${TMPDIR}/liq_Het_Coincidence.res" ] && cp ${TMPDIR}/liq_Het_Coincidence.res ${SCRDIR}
# [ -d "${TMPDIR}/poly_Het_Coincidence.res" ] && cp ${TMPDIR}/poly_Het_Coincidence.res ${SCRDIR}
# [ -d "${TMPDIR}/liq_PRE_Coincidence.res" ] && cp ${TMPDIR}/liq_PRE_Coincidence.res ${SCRDIR}
# [ -d "${TMPDIR}/poly_PRE_Coincidence.res" ] && cp ${TMPDIR}/poly_PRE_Coincidence.res ${SCRDIR}
# [ -d "${TMPDIR}/polyAniso.res" ] && cp ${TMPDIR}/polyAniso.res ${SCRDIR}
# [ -d "${TMPDIR}/polyGyration.res" ] && cp ${TMPDIR}/polyGyration.res ${SCRDIR}
# [ -d "${TMPDIR}/polyHetMSD.res" ] && cp ${TMPDIR}/polyHetMSD.res ${SCRDIR}
# [ -d "${TMPDIR}/polyHomMSD.res" ] && cp ${TMPDIR}/polyHomMSD.res ${SCRDIR}
# [ -d "${TMPDIR}/polyPREMSD.res" ] && cp ${TMPDIR}/polyPREMSD.res ${SCRDIR}
# [ -d "${TMPDIR}/liqMean.res" ] && cp ${TMPDIR}/liqMean.res ${SCRDIR}
# [ -d "${TMPDIR}/liqSTD.res" ] && cp ${TMPDIR}/liqSTD.res ${SCRDIR}

cp ${TMPDIR}/liqDropNum.res ${SCRDIR}
cp ${TMPDIR}/polyHetMSD.res ${SCRDIR}
cp ${TMPDIR}/polyHomMSD.res ${SCRDIR}
cp ${TMPDIR}/polyPREMSD.res ${SCRDIR}

# Archive output files to home directory
tar -czf ${DATDIR}/${STEP}.tar.gz ${TMPDIR}/

# Clean scratch
rm -rf ${TMPDIR}
