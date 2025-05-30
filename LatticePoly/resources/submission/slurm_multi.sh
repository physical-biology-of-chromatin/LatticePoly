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
sed -e "${DIRSUB}" < /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/submission/${TEMP}/input.cfg > ${TMPDIR}/input.cfg


# Run
./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out

#Analysis


python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  255 5 20
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  100 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  200 2100

#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  300 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  400 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  550 1 50
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  550 2 50
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  550 3 50
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  550 4 50
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  550 6 50
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  500 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  700 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}   2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  520 5 80
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  600 4 100

#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  600 4 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  700 4 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  800 4 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  900 4 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  1000 4 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  1200 4 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  1400 4 100
#python3 /home/ddasaro//Yeast_full_genome/LatticePoly/LatticePoly/resources/Forks_n_percentage_full_genome.py ${TMPDIR} 500
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  700 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  700 6 100

#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 4 400
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  750 4 20


#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  520 4 80
#MonomerDist_HiC_full_genome_cis.py
python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_trans.py ${TMPDIR}  550 4 50
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_trans.py ${TMPDIR}  550 3 50
python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_cis.py ${TMPDIR}  550 4 50
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_cis.py ${TMPDIR}  550 3 50


#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  900 3 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  1100 3 100

#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 4 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 5 100

#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 3 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 4 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 5 100


python3 /home/ddasaro//Yeast_full_genome/LatticePoly/LatticePoly/resources/Forks_n_percentage_full_genome.py ${TMPDIR} 0
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 0
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 1
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 2
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 3
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 4
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 5
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 6
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 7
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 8
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 9
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 10
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 11
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 12
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 13
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 14
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 15
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 16
python3 /home/ddasaro/LatticePoly/LatticePoly/resources/Forksnumber_chrom.py ${TMPDIR} 0 17

# Move SLURM output files to data directory
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${TMPDIR}
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${TMPDIR}

# Create data directory on local disk
#[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}
#mkdir /Xnfs/physbiochrom/ddasaro/year3/October24/full_genome_S/
# Archive output files to home directory
tar --transform "s|^|${OUTDIR}/|" -czf /Xnfs/physbiochrom/ddasaro/year4/simulations_paper_cohesion/${OUTDIR}.tar.gz -C ${TMPDIR} .
#tar -xzf ${DATDIR}/${OUTDIR}.tar.gz -C /Xnfs/lbmcdb/Jost_team/ddasaro/year2/April23/Ring/100
#rm -rf {DATDIR}/${OUTDIR}.tar.gz

# Clean scratch
#rm -rf ${TMPDIR}


