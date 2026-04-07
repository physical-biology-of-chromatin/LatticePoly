import os
import sys
import numpy as np
import time
import subprocess

if __name__ == "__main__":
        if len(sys.argv) != 2:
                print("\033[1;31mUsage is %s outputDir \033[0m" % sys.argv[0])
                #sys.exit()
                outputDir = "Sim1"
        else:
                outputDir = sys.argv[1]

SCRIPTDIR = os.path.join(os.getcwd(),'resources/submission')

def string_to_list(input):
        if not ',' in input:
                return([input.strip()])
        else:
                return([i.strip() for i in input.split(",")])

dict_parameters = {}

file = open("resources/submission/input_slurm.cfg",'r')
for line in file.readlines():
        if line.split(' = ')[0] == "Nstat":
                Nstat = int(line.split(' = ')[1])
        else:
                dict_parameters[line.split(' = ')[0]] = string_to_list(line.split(' = ')[1])
file.close()


nb_job = int(np.prod(np.array([len(values) for keys,values in dict_parameters.items()])))

nb_batch = nb_job//48
if nb_job%48 != 0:
        raise IOError("Number of job (%d) is not a multiple of 48" % nb_job)


dict_mapping_parameters = {}

for keys,values in dict_parameters.items():
        if len(values) != 1:
                dict_mapping_parameters[keys] = len(values)

dict_parameters_for_job = {keys : [] for keys in dict_mapping_parameters.keys()}

compteur = 1

for keys, values in dict_mapping_parameters.items():
        for job in range(nb_job):
                dict_parameters_for_job[keys].append(dict_parameters[keys][(job//compteur)%values])
        compteur *= values
batch = 0

for n in range(Nstat):
        for job in range(nb_job):

                file = open("resources/submission/slurm_sweep_tmp.sh", "w")
                
                file.write("#!/bin/bash\n##\n##  slurm_sweep_tmp.sh\n##  LatticePoly\n")
                file.write(f"##\n##  Created by ppuel on {time.localtime().tm_mday}/{time.localtime().tm_mon}/{time.localtime().tm_year}\n")
                file.write("##  Copyright © 2022 ENS Lyon. All rights reserved.\n##\n\n")
                
                file.write("#SBATCH -o %A_%a.out\n#SBATCH -e %A_%a.err\n\n")

                file.write("# Output directory\n")
                file.write(f"OUTPUTDIR={outputDir}\n")

                file.write("\n# Temporary directory\n")
                file.write("TEMPORARYDIR=/home/${LOGNAME}/Documents/Simulation/tmp/tmp\n")

                file.write("\n# Associated scratch directory\n")
                file.write("SCRATCHDIR=/home/${LOGNAME}/Documents/Simulation/tmp/data\n")

                file.write("\n# Script (relative) path\n")
                file.write(f"SCRIPTDIR={SCRIPTDIR}\n")
                
                file.write("\n# Data directory\n")
                file.write("XNFSDIR=/home/${LOGNAME}/Documents/Simulation/tmp/archive\n")

                file.write("\n# Relative path to code root directory\n")
                file.write("ROOTDIR=${SCRIPTDIR}/../..\n\n")
                
                file.write("# Set working directory to root\ncd ${ROOTDIR}\n\n")
                file.write("LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTDIR/lib\n")
                file.write("PYTHONPATH=$PYTHONPATH:$ROOTDIR/resources/h5py\n")
                file.write("export LD_LIBRARY_PATH\nexport PYTHONPATH\n\n")
                
                file.write("# Executable path\nEXEC=bin/lat\n\n")
                
                file.write("#Values of the parameters\n")

                file.write(f"SLURM_ARRAY_TASK_ID={job}\n")

                for keys, values in dict_parameters.items():
                        if len(values) == 1:
                                file.write(f"{keys.upper()}={values[0]}\n")
                        else:
                                file.write(f"\n{keys.upper()}_ARRAY=({" ".join(dict_parameters_for_job[keys][batch*48:(batch+1)*48])})\n")
                                file.write(f"{keys.upper()}=$"+"{"+f"{keys.upper()}_ARRAY[$"+"[$SLURM_ARRAY_TASK_ID-1]]}\n\n")

                file.write("# Temporary directory on tmp\n")

                temporary_path = "${TEMPORARYDIR}/${OUTPUTDIR}"
                for keys, values in dict_mapping_parameters.items():
                        temporary_path += f"_{keys.upper()}_$"+"{"+f"{keys.upper()}"+"}"
                temporary_path += f"_N_{n}"
        
                file.write(f"TMPDIR={temporary_path}\n\n")
                
                file.write("\n# H5 File path\n")
                file.write("FILEPATH=${TMPDIR}/traj.h5\n")

                file.write("# Create temporary directory if necessary\n[ ! -d \"${TMPDIR}\" ] && mkdir -p ${TMPDIR}\n\n")

                file.write("# Scratch directory on scratch\n")

                scratch_path = "${SCRATCHDIR}/${OUTPUTDIR}"
                for keys, values in dict_mapping_parameters.items():
                        scratch_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
                scratch_path += f"/N/{n}"
        
                file.write(f"SCRDIR={scratch_path}\n\n")
                
                file.write("# Create scratch directory if necessary\n[ ! -d \"${SCRDIR}\" ] && mkdir -p ${SCRDIR}\n\n")
                

                file.write("# Data directory on Xnfs\n")

                data_path = "${XNFSDIR}/${OUTPUTDIR}"
                for keys, values in dict_mapping_parameters.items():
                        data_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
                data_path += f"/N/{n}"

                file.write(f"DATDIR={data_path}\n\n")
                
                file.write("# Create data directory if necessary\n[ ! -d \"${DATDIR}\" ] && mkdir -p ${DATDIR}\n\n")
               

                file.write("# Substitution strings\n")
                
                file.write("DIRSUB=\"s|\\(outputDir[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${TMPDIR} ;|;\"\n")
                file.write("FILSUB=\"s|\\(H5filePath[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${FILEPATH} ;|;\"\n")
                
                for keys in dict_parameters.keys():
                        file.write(f"{keys.upper()}SUB=\"s|\\("+f"{keys}"+"[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${"+f"{keys.upper()}"+"} ;|;\"\n")
                
                
                file.write("\n# Copy input configuration file to output directory, substituting paths and parameter values\n")
                
                sed_string = "sed -e \"${DIRSUB}\"\"${FILSUB}\""
                for keys in dict_parameters.keys():
                        sed_string += "\"${"+f"{keys.upper()}SUB"+"}\""
                
                file.write(sed_string+" < data/input.cfg > ${TMPDIR}/input.cfg\n")

                file.write("\n# Run\n./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out\n")
                
                file.write("\n# Perform post-processing analyses\n")
                file.write(".venv/bin/python3 resources/h5py/LiqDensity.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out\n")
                file.write(".venv/bin/python3 resources/h5py/LiqCluster.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out\n")
                file.write(".venv/bin/python3 resources/h5py/LiqMSD.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out\n")
                file.write(".venv/bin/python3 resources/h5py/LiqPolyCoincidence.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out\n")
                file.write(".venv/bin/python3 resources/h5py/PolyGyration.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out\n")
                file.write(".venv/bin/python3 resources/h5py/PolyMSD.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out\n\n")
                
                file.write("\n# Move processed output files to scratch directory\n")
                file.write("cp ${TMPDIR}/input.cfg ${SCRDIR}\n")
                file.write("cp ${TMPDIR}/log.out ${SCRDIR}\n")
                file.write("cp ${TMPDIR}/process.out ${SCRDIR}\n")
                file.write("cp ${TMPDIR}/process.h5 ${SCRDIR}\n")
                
                file.write("\n# Move all files to XNFS directory\n")
                file.write("mv ${FILEPATH} ${DATDIR}/\n")
                file.write("mv ${TMPDIR}/input.cfg ${DATDIR}\n")
                file.write("mv ${TMPDIR}/log.out ${DATDIR}\n")
                file.write("mv ${TMPDIR}/process.out ${DATDIR}\n")
                file.write("mv ${TMPDIR}/process.h5 ${DATDIR}\n")

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
                JOBNAME = ""
                for keys in dict_mapping_parameters.keys():
                        JOBNAME += (keys+"_") 
                JOBNAME += str(nb_batch*n+batch)

                # sbatch command
                command = f"bash /{SCRIPTDIR}/slurm_sweep_tmp.sh"              
                print(command)
                subprocess.run(command, shell=True, executable="/bin/bash")
                print("why")
                


 
