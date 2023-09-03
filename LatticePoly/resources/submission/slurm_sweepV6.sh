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

ID_FIELD=${SLURM_ARRAY_TASK_ID}

FIELD=data/spinField_${NUMFIELD}.in

J_STEP=$((${NUM_STEP}-1))

EVARRAY=("0.0" "0.5" "1.0" "1.5" "2.0" "2.5" "3.0" "3.5" "4.0")

JLPARRAY=("0.0" "0.05" "0.1" "0.15" "0.2" "0.25" "0.3" "0.35")

# jlp_index = i%16//2 
# ev_index = j*3+i//16
# l = i%2

JLP=${JLPARRAY[$((((ID_FIELD-1)%16)/2))]}
# echo ${JLP} >> ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out

EV=${EVARRAY[$((J_STEP*3+(ID_FIELD-1)/16))]}
# echo ${EV} >> ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out

STEP=$(((ID_FIELD-1)%2))


# Temporary directory on tmp
TMPDIR=spinField_${NUMFIELD}_JLP_${JLP}_EV_${EV}_${STEP}

TMPDIR=${TEMPORARYDIR}/${OUTPUTDIR}_${GENE}_${TMPDIR}

# Create temporary directory if necessary
[ ! -d "${TMPDIR}" ] && mkdir -p ${TMPDIR}

# Scratch directory on scratch
SCRDIR=spinField_${NUMFIELD}/JLP/${JLP}/EV/${EV}/${STEP}

SCRDIR=${SCRATCHDIR}/${OUTPUTDIR}/${GENE}/${SCRDIR}

# Create scratch directory if necessary
[ ! -d "${SCRDIR}" ] && mkdir -p ${SCRDIR}

# Data directory on Xnfs
DATDIR=spinField_${NUMFIELD}/JLP/${JLP}/EV/${EV}

DATDIR=${XNFSDIR}/${OUTPUTDIR}/${GENE}/${DATDIR}

# Create data directory if necessary
[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}



# Substitution strings
DIRSUB="s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${TMPDIR} ;|;"
JLPSUB="s|\(Jlp[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JLP} ;|;"

EVSUB="s|\(EV[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${EV} ;|;"
FIELDSUB="s|\(fieldPath[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${FIELD} ;|;"


# Copy input configuration file to output directory, substituting paths and parameter values
sed -e "${DIRSUB}""${JLPSUB}""${EVSUB}""${FIELDSUB}" < data/input.cfg > ${TMPDIR}/input.cfg


# Run
./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out

# Perform post-processing analyses
python3 resources/DistanceMap.py ${TMPDIR} -1 10 >> ${TMPDIR}/process.out
python3 resources/PolyGyration.py ${TMPDIR} -1 >> ${TMPDIR}/process.out
python3 resources/PolyMSD.py ${TMPDIR} -1 >> ${TMPDIR}/process.out
python3 resources/Probe_distance.py ${TMPDIR} -1 ${GENE} >> ${TMPDIR}/process.out



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

cp ${TMPDIR}/polyAniso.res ${SCRDIR}
cp ${TMPDIR}/polyGyration.res ${SCRDIR}
cp ${TMPDIR}/polyHetMSD.res ${SCRDIR}
cp ${TMPDIR}/polyHomMSD.res ${SCRDIR}
cp ${TMPDIR}/polyPREMSD.res ${SCRDIR}
cp ${TMPDIR}/Probe_distance.res ${SCRDIR}

# Archive output files to home directory
tar -czf ${DATDIR}/${STEP}.tar.gz ${TMPDIR}/

# Clean scratch
rm -rf ${TMPDIR}
