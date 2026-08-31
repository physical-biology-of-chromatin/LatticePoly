##
##  submit_slurm_Trajectory.py
##  LatticePoly
##
##  Created by ppuel on 29/10/2024.
##  Copyright © 2024 ENS Lyon. All rights reserved.
##

import os, sys, time, subprocess
import numpy as np

if __name__ == "__main__":
        if len(sys.argv) != 2:
                print("\033[1;31mUsage is %s outputDir \033[0m" % sys.argv[0])
                sys.exit()
        else:
                outputDir = sys.argv[1]
                outputDir_format = outputDir.replace("/","_")


def string_to_list(input):
        if not ',' in input:
                return([input.strip()])
        else:
                return([i.strip() for i in input.split(",")])


SCRIPTDIR = os.path.join(os.getcwd(),'resources/submission')

dict_parameters = {}

is_poly = True

file = open("resources/submission/input_slurm.cfg",'r')
for line in file.readlines():
        if line.split(' = ')[0] == "Nstat":
                dict_parameters["N"] = [str(i+40) for i in range(int(line.split(' = ')[1]))]
        else:
                dict_parameters[line.split(' = ')[0]] = string_to_list(line.split(' = ')[1])
                if line.split(' = ')[1].strip() == 'data/toy_domain.in':
                        is_poly = False

file.close()

os.makedirs(os.path.join("/Xnfs/physbiochrom/ppuel/data/",outputDir), exist_ok=True)

print(subprocess.run(f"cp /home/ppuel/Simulation/LatticePoly/LatticePoly/resources/submission/input_slurm.cfg {os.path.join('/Xnfs/physbiochrom/ppuel/data/',outputDir)}", shell=True, executable="/bin/bash"))

nb_job = int(np.prod(np.array([len(values) for values in dict_parameters.values()])))

nb_batch = nb_job//48
if nb_job%48 != 0:
        raise IOError("Number of job (%d) is not a multiple of 48" % nb_job)


dict_mapping_parameters = {}

for keys,values in dict_parameters.items():
        if len(values) != 1:
                dict_mapping_parameters[keys] = len(values)

dict_parameters_for_job = {keys : [] for keys in dict_mapping_parameters.keys()}
dict_parameters_for_job['N'] = []

compteur = 1

for keys, values in dict_mapping_parameters.items():
        for job in range(nb_job):
                dict_parameters_for_job[keys].append(dict_parameters[keys][(job//compteur)%values])
        compteur *= values

compteur_batch = 0

for batch in range(nb_batch):

        compteur_batch += 1

        file = open(f"resources/submission/tmp/slurm_sweep_{outputDir_format}_{batch}_{compteur_batch}.sh", "w")

        file.write("#!/bin/bash\n##\n##  slurm_sweep_tmp.sh\n##  LatticePoly\n")
        file.write(f"##\n##  Created by ppuel on {time.localtime().tm_mday}/{time.localtime().tm_mon}/{time.localtime().tm_year}\n")
        file.write(f"##  Copyright © {time.localtime().tm_year} ENS Lyon. All rights reserved.\n##\n\n")

        file.write("#SBATCH -o %A_%a.out\n#SBATCH -e %A_%a.err\n\n")

        file.write("# Output directory\n")
        file.write(f"OUTPUTDIR={outputDir}\n")

        file.write("\n# Temporary directory\n")
        file.write("TEMPORARYDIR=/tmp/${LOGNAME}\n")

        file.write("\n# Associated scratch directory\n")
        file.write("SCRATCHDIR=/scratch/Cascade/${LOGNAME}/data\n")

        file.write("\n# Script (relative) path\n")
        file.write(f"SCRIPTDIR={SCRIPTDIR}\n")

        file.write("\n# Data directory\n")
        file.write("XNFSDIR=/Xnfs/physbiochrom/${LOGNAME}/data\n")

        file.write("\n# Relative path to code root directory\n")
        file.write("ROOTDIR=${SCRIPTDIR}/../..\n\n")

        file.write("# Set working directory to root\ncd ${ROOTDIR}\n\n")
        file.write("LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTDIR/lib\n")
        file.write("PYTHONPATH=$PYTHONPATH:$ROOTDIR/resources/h5py\n")
        file.write("export LD_LIBRARY_PATH\nexport PYTHONPATH\n\n")

        file.write("# Executable path\nEXEC=bin/lat\n\n")

        file.write("# Values of the parameters\n")

        for keys, values in dict_parameters.items():
                if len(values) == 1:
                        file.write(f"{keys.upper()}={values[0]}\n")
                else:
                        file.write(f"\n{keys.upper()}_ARRAY=({' '.join(dict_parameters_for_job[keys][batch*48:(batch+1)*48])})\n")
                        file.write(f"{keys.upper()}=$"+"{"+f"{keys.upper()}_ARRAY[$"+"[$SLURM_ARRAY_TASK_ID-1]]}\n\n")
        
        file.write("# Temporary directory on tmp\n")

        temporary_path = "${TEMPORARYDIR}/${OUTPUTDIR}"
        for keys, values in dict_mapping_parameters.items():
                if keys != 'N':
                        temporary_path += f"_{keys.upper()}_$"+"{"+f"{keys.upper()}"+"}"
        temporary_path += "_N_${N}"

        file.write(f"TMPDIR={temporary_path}\n\n")

        file.write("\n# H5 File path\n")
        file.write("FILEPATH=${TMPDIR}/traj.h5\n")

        file.write("# Create temporary directory if necessary\n[ ! -d \"${TMPDIR}\" ] && mkdir -p ${TMPDIR}\n\n")

        file.write("# Scratch directory on scratch\n")

        scratch_path = "${SCRATCHDIR}/${OUTPUTDIR}"
        for keys, values in dict_mapping_parameters.items():
                if keys != "N":
                        scratch_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
        scratch_path += "/N/${N}"

        file.write(f"SCRDIR={scratch_path}\n\n")

        file.write("# Create scratch directory if necessary\n[ ! -d \"${SCRDIR}\" ] && mkdir -p ${SCRDIR}\n\n")


        file.write("# Data directory on Xnfs\n")

        data_path = "${XNFSDIR}/${OUTPUTDIR}"
        for keys, values in dict_mapping_parameters.items():
                if keys != "N":
                        data_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
        data_path += "/N/${N}"

        file.write(f"DATDIR={data_path}\n\n")

        file.write("# Create data directory if necessary\n[ ! -d \"${DATDIR}\" ] && mkdir -p ${DATDIR}\n\n")


        file.write("# Substitution strings\n")

        file.write("DIRSUB=\"s|\\(outputDir[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${TMPDIR} ;|;\"\n")
        file.write("FILSUB=\"s|\\(H5filePath[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${FILEPATH} ;|;\"\n")

        for keys in dict_parameters.keys():
                if not(keys in ["mode", "exponant"]):
                        file.write(f"{keys.upper()}SUB=\"s|\\("+f"{keys}"+"[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${"+f"{keys.upper()}"+"} ;|;\"\n")


        file.write("\n# Copy input configuration file to output directory, substituting paths and parameter values\n")

        sed_string = "sed -e \"${DIRSUB}\"\"${FILSUB}\""
        for keys in dict_parameters.keys():
                if not(keys in ["mode", "exponant"]):
                        sed_string += "\"${"+f"{keys.upper()}SUB"+"}\""

        file.write(sed_string+" < data/input.cfg > ${TMPDIR}/input.cfg\n")

        file.write("\n# Run\n./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out\n")

        file.write("\n# Perform post-processing analyses\n")
        file.write(".venv/bin/python3 resources/h5py/LiqDensity.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
        file.write(".venv/bin/python3 resources/h5py/LiqCluster.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
        file.write(".venv/bin/python3 resources/h5py/LiqMSD.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
        file.write(".venv/bin/python3 resources/h5py/LiqDroplet.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
        
        if is_poly:
                file.write(".venv/bin/python3 resources/h5py/PolyMSD.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
                file.write(".venv/bin/python3 resources/h5py/PolyGyration.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
                file.write(".venv/bin/python3 resources/h5py/LiqPolyCoM.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
                file.write(".venv/bin/python3 resources/h5py/PolyLiqDensity.py ${TMPDIR} 17 >> ${TMPDIR}/process.out\n")
                file.write(".venv/bin/python3 resources/h5py/PolyTad.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
                
        file.write("\n# Move processed output files to scratch directory\n")
        file.write("cp ${FILEPATH} ${SCRDIR}\n")
        file.write("cp ${TMPDIR}/process.h5 ${SCRDIR}\n")

        file.write("\n# Move all files to XNFS directory\n")
        file.write("mv ${FILEPATH} ${DATDIR}/\n")
        file.write("mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${DATDIR}\n")
        file.write("mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${DATDIR}\n")
        file.write("mv ${TMPDIR}/input.cfg ${DATDIR}\n")
        file.write("mv ${TMPDIR}/log.out ${DATDIR}\n")
        file.write("mv ${TMPDIR}/process.out ${DATDIR}\n")
        file.write("mv ${TMPDIR}/process.h5 ${DATDIR}\n")
        # file.write("mv ${TMPDIR}/liq_graph ${DATDIR}\n")
        file.write("mv ${TMPDIR}/liq_droplets.pickle ${DATDIR}\n")
        file.write("mv ${TMPDIR}/liq_simple_droplets.pickle ${DATDIR}\n")

        file.write("\n# Clean scratch\n")
        file.write("rm -rf ${TMPDIR}\n")

        file.close()

        # Max. walltime
        WTIME = "8-00:00:00"

        # Partition
        QUEUE = "Cascade"

        # Max. memory per task
        MAXMEM = "1G"

        # sbatch arguments
        QARGS=f"--ntasks=1 --mem={MAXMEM} -a 1-48 -t {WTIME} -p {QUEUE}"

        # Job Name
        JOBNAME = outputDir_format+str(batch)

        # sbatch command
        command = f"sbatch {QARGS} -J {JOBNAME} {SCRIPTDIR}/tmp/slurm_sweep_{outputDir_format}_{batch}_{compteur_batch}.sh"
        subprocess.run(command, shell=True, executable="/bin/bash")
        