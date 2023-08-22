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





python /home/ddasaro//LatticePoly/LatticePoly/resources/Poly_Rcmdiff_SCs.py ${TMPDIR} 2000
python /home/ddasaro//LatticePoly/LatticePoly/resources/Gyration_moments.py ${TMPDIR} 2000
python /home/ddasaro//LatticePoly/LatticePoly/resources/PolyGyration_repl_chr1.py ${TMPDIR} 2000









#python3 /home/ddasaro//LatticePoly/LatticePoly/resources/PolyGyration_repl_chr1.py  ${TMPDIR}  25000
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 10  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 20  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 40  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 80  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   200 500 100  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 200  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 300  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 400  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 500
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 600  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 700  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 800  
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 900 
python3 /home/ddasaro//LatticePoly/LatticePoly/resources/Distance_mon.py ${TMPDIR}   2000 500 999




python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_trans_intervals_cool.py ${TMPDIR} 20000 0 2 
python /home/ddasaro//LatticePoly/LatticePoly/resources/MonomerDist_HiC_M_cis_intervals_cool.py ${TMPDIR} 20000 0 2

python /home/ddasaro//LatticePoly/LatticePoly/resources/ReplicationAnalysis.py ${TMPDIR} 2000


# Move SLURM output files to data directory
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${TMPDIR}
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${TMPDIR}

# Create data directory on local disk
#[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Archive output files to home directory
tar --transform "s|^|${OUTDIR}/|" -czf /Xnfs/lbmcdb/Jost_team/ddasaro/year2/July23/bubble_long/null/${OUTDIR}.tar.gz -C ${TMPDIR} .
#tar -xzf ${DATDIR}/${OUTDIR}.tar.gz -C /Xnfs/lbmcdb/Jost_team/ddasaro/year2/April23/Ring/100
#rm -rf {DATDIR}/${OUTDIR}.tar.gz

# Clean scratch
#rm -rf ${TMPDIR}


