#$ -S /bin/bash

### Max. walltime
#$ -l h_rt=96:00:00

### Queue name(s)
#$ -q CLG6242deb384A,CLG6242deb384C,CLG5218deb192A,CLG5218deb192B,CLG5218deb192C,CLG5218deb192D

### Load the user environment for SGE
#$ -cwd

### Export environment variables to all runtime nodes
#$ -V

# Run N parallel jobs, labelled from 1 to N
#$ -t 1-${N}


# Work from code root directory
ROOTDIR=${HOMEDIR}/LatticePoly/LatticePoly

cd ${ROOTDIR}

# Executable path
EXEC=bin/lat

# Scratch directory
SCRATCHDIR=/scratch/Bio/${LOGNAME}

# Set parameter value by linear interpolation based on job ID
VAL=`echo ${MINVAL} ${MAXVAL} ${N} ${SGE_TASK_ID} | awk '{printf("%.3f\n", $1+($2-$1)/($3-1)*($4-1))}'`

# Ouput directory
OUTDIR=${SCRATCHDIR}/LatticeData/${PARAM}_${VAL}

# Create directory if necessary
if [[ ! -d "${OUTDIR}" ]]
then
	mkdir -p ${OUTDIR}
fi

# Copy input configuration file to output directory, substituting paths and parameter values
sed -e 's/\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)/\1${OUTDIR} ;/;s/\(${PARAM}[[:space:]]*=[[:space:]]*\)\(.*;\)/\1${VAL} ;/' < data/input.cfg > ${OUTDIR}/.input.cfg

# Run
./${EXEC} ${OUTDIR}/.input.cfg > ${OUTDIR}/log.out

# Copy output files to home directory
cp -r ${OUTDIR} data

# Clear scratch
rm -r ${OUTDIR}
