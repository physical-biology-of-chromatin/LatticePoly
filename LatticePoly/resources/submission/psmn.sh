#$ -S /bin/bash

### Max. walltime
#$ -l h_rt=96:00:00
### Queue name(s)
#$ -q CLG6242deb384A,CLG6242deb384C,CLG5218deb192A,CLG5218deb192B,CLG5218deb192C,CLG5218deb192D
### Load the user environment for SGE
#$ -cwd
### Export environment variables to all runtime nodes
#$ -V


# Code root directory
ROOTDIR=${HOME}/LatticePoly/LatticePoly

# Set working directory to root
cd ${ROOTDIR}

# Set parameter values by linear interpolation between MIN_VAL and MAX_VAL based on task ID
VAL=$(echo ${MIN_VAL} ${MAX_VAL} ${SGE_TASK_FIRST} ${SGE_TASK_LAST} ${SGE_TASK_ID} | awk '{printf("%.3f\n", $1+($2-$1)/($4-$3)*($5-$3))}')

# Executable path
EXEC=bin/lat

# Scratch directory
SCRATCHDIR=/scratch/Bio/${LOGNAME}

# Ouput directory
OUTDIR=${SCRATCHDIR}/LatticeData/${PARAM}_${VAL}

# Create output directory if necessary
if [ ! -d "${OUTDIR}" ]; then
	mkdir -p ${OUTDIR}
fi

# Copy input configuration file to output directory, substituting directory paths and parameter values
sed -e "s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${OUTDIR} ;|;s|\(${PARAM}[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${VAL} ;|" < data/input.cfg > ${OUTDIR}/.input.cfg

# Run
./${EXEC} ${OUTDIR}/.input.cfg > ${OUTDIR}/log.out

# Move SGE output files to data directory
mv ${JOB_NAME}.e${JOB_ID}.${SGE_TASK_ID} ${OUTDIR}
mv ${JOB_NAME}.o${JOB_ID}.${SGE_TASK_ID} ${OUTDIR}

# Create data folder in home directory
if [ ! -d "data/${PARAM}" ]; then
	mkdir -p data/${PARAM}
fi

# Archive output files to home directory
tar zcvf data/${PARAM}/${VAL}.tar.gz ${OUTDIR}

# Clean scratch
rm -r ${OUTDIR}
