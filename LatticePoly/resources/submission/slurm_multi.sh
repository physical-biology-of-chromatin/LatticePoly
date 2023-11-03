#!/bin/bash

##
##  slurm.sh
##  LatticePoly
##
##  Created by ddasaro on the model of mtortora on 20/06/2022
##  Copyright © 2022 ENS Lyon. All rights reserved.
##


#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# Relative path to code root directory
ROOTDIR=${SCRIPTDIR}/../..
echo "The present scriptdir directory is ${SCRIPTDIR}";



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
[ ! -d "${SCRATCHDIR}/${LOGNAME}/LatticeData/${SLURM_JOB_NAME}/" ] && mkdir -p ${SCRATCHDIR}/${LOGNAME}/LatticeData/${SLURM_JOB_NAME}/
mkdir ${SCRATCHDIR}/${LOGNAME}/LatticeData/${SLURM_JOB_NAME}/${OUTDIR}
TMPDIR=${SCRATCHDIR}/${LOGNAME}/LatticeData/${SLURM_JOB_NAME}/${OUTDIR}

# Create output directory if necessary
#[ ! -d "${TMPDIR}" ] && mkdir -p ${TMPDIR}

# Substitution strings
DIRSUB="s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${TMPDIR} ;|;"

# Copy input configuration file to output directory, substituting paths and parameter values
sed -e "${DIRSUB}" < /home/ddasaro/LatticePoly/LatticePoly/resources/submission/${TEMP}/input.cfg > ${TMPDIR}/input.cfg


# Run
./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out

#Analysis








python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_S_cool.py ${TMPDIR} 500 1 100 
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_S_cool.py ${TMPDIR} 500 2 100 

python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_S_cool.py ${TMPDIR} 900 1 100 
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_S_cool.py ${TMPDIR} 900 2 100 

python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_S_cool.py ${TMPDIR} 1300 1 100 
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_S_cool.py ${TMPDIR} 1300 2 100 

python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_S_cool.py ${TMPDIR} 1700 1 100 
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_S_cool.py ${TMPDIR} 1700 2 100 








python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Forksnumber.py ${TMPDIR} 500
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/PolyMSDchromatid1.py  ${TMPDIR}  5010
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py ${TMPDIR}   5010 2
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py ${TMPDIR}   5010 3 
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py ${TMPDIR}   5010 4 

#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_cool.py ${TMPDIR}   500 2 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_cool.py ${TMPDIR}   500 4 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_cool.py ${TMPDIR}  500 2 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_cool.py ${TMPDIR}   500 4 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_cool.py ${TMPDIR}   500 3 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_cool.py ${TMPDIR}  500 3 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_cool.py ${TMPDIR}   500 5 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_cool.py ${TMPDIR}   500 5 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_all_cool.py ${TMPDIR}   500 2 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_all_cool.py ${TMPDIR}   500 5 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_all_cool.py ${TMPDIR}   500 3 chrIV
#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_all_cool.py ${TMPDIR}   500 4 chrIV
# Move SLURM output files to data directory
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${TMPDIR}
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${TMPDIR}

# Create data directory on local disk
#[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Archive output files to home directory
tar --transform "s|^|${OUTDIR}/|" -czf /Xnfs/lbmcdb/Jost_team/ddasaro/year3/October23/cohesion_S/${OUTDIR}.tar.gz -C ${TMPDIR} .
#tar -xzf ${DATDIR}/${OUTDIR}.tar.gz -C /Xnfs/lbmcdb/Jost_team/ddasaro/year2/April23/Ring/100
#rm -rf {DATDIR}/${OUTDIR}.tar.gz

# Clean scratch
#rm -rf ${TMPDIR}



