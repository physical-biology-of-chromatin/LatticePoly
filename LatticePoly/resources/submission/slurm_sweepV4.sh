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

#QVARS="OUTPUTDIR=$1,PARAM=$2,MIN_VAL=$3,MAX_VAL=$4,SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR}"
#QVARS2="PARAM2=$5,MIN_VAL2=$6,MAX_VAL2=$7,STEP=$s,ITER=$i"  {printf("%.4f\n", $1+($2-$1)/($4/4-$3)*($5+($6-1)*4)}

TASK_ID_RESCALLED=$(((SLURM_ARRAY_TASK_ID-1)/12))

VAL=$(echo ${MIN_VAL} ${MAX_VAL} ${SLURM_ARRAY_TASK_MIN} ${SLURM_ARRAY_TASK_MAX} ${SLURM_ARRAY_TASK_ID} | awk '{printf("%.4f\n", $1+($2-$1)/($4/4-$3)*(($5-$3)%12))}')
VAL2=$(echo ${MIN_VAL2} ${MAX_VAL2} ${SLURM_ARRAY_TASK_MIN} ${SLURM_ARRAY_TASK_MAX} ${TASK_ID_RESCALLED} ${STEP} | awk '{printf("%.4f\n", $1+($2-$1)/($4/4-$3)*($5+($6-1)*4))}')

# Temporary directory on tmp
TMPDIR=${PARAM}_${VAL}_${ITERATION}
[ ! -z "${PARAM2}" ] && TMPDIR=${PARAM2}_${VAL2}_${TMPDIR}

TMPDIR=${TEMPORARYDIR}/${OUTPUTDIR}_${TMPDIR}

# Create temporary directory if necessary
[ ! -d "${TMPDIR}" ] && mkdir -p ${TMPDIR}

# Scratch directory on scratch
SCRDIR=${PARAM}/${VAL}/${ITERATION}
[ ! -z "${PARAM2}" ] && SCRDIR=${PARAM2}/${VAL2}/${SCRDIR}

SCRDIR=${SCRATCHDIR}/${OUTPUTDIR}/${SCRDIR}

# Create scratch directory if necessary
[ ! -d "${SCRDIR}" ] && mkdir -p ${SCRDIR}

# Data directory on Xnfs
DATDIR=${PARAM}/${VAL}
[ ! -z "${PARAM2}" ] && DATDIR=${PARAM2}/${VAL2}/${DATDIR}

DATDIR=${XNFSDIR}/${OUTPUTDIR}/${DATDIR}

# Create data directory if necessary
[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Substitution strings
DIRSUB="s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${TMPDIR} ;|;"
VALSUB="s|\(${PARAM}[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${VAL} ;|;"

[ ! -z "${PARAM2}" ] && VAL2SUB="s|\(${PARAM2}[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${VAL2} ;|;"

# Copy input configuration file to output directory, substituting paths and parameter values
sed -e "${DIRSUB}""${VALSUB}""${VAL2SUB}" < data/input.cfg > ${TMPDIR}/input.cfg

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
tar -czf ${DATDIR}/${ITERATION}.tar.gz ${TMPDIR}/

# Clean scratch
rm -rf ${TMPDIR}
