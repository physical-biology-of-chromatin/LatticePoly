##
##  submit_slurm_process.py
##  LatticePoly
##
##  Created by ppuel on 29/10/2024.
##  Copyright © 2024 ENS Lyon. All rights reserved.
##

import os
import sys
import numpy as np
import time
import subprocess

if __name__ == "__main__":
        if len(sys.argv) < 3:
                print("\033[1;31mUsage is %s experience program1 data1 program2 data2 ...\033[0m" % sys.argv[0])
                sys.exit()
        else:
                experience = sys.argv[1]

                dict_program = {}

                arg = 2
                tmp_arg = None
                while arg<len(sys.argv):
                        if "resources/" in sys.argv[arg]:
                                dict_program[sys.argv[arg]] = []
                                tmp_arg = arg
                                arg+=1
                        else:
                                dict_program[sys.argv[tmp_arg]].append(sys.argv[arg])
                                arg+=1
                
                rm_last_process = 0
        
        XnfsDir = "/Xnfs/physbiochrom/ppuel/data/"
        ScratchDir = "/scratch/Cascade/ppuel/data/"

        verif = 0

        for tmp_experience in os.listdir(XnfsDir):
                if tmp_experience.split("_")[0] == f'EXP{experience}':
                        outputDir = tmp_experience
                        verif += 1

        if verif != 1:
                print(f"The experience {experience} is not found or found multiple time")
                sys.exit()      
        
        
        outputDir_format = outputDir.replace("/","_")


def string_to_list(input):
        if not ',' in input:
                return([input.strip()])
        else:
                return([i.strip() for i in input.split(",")])
        
def arg_to_variable_slurm(list_arg, dict_parameters):
        string = ""

        for arg in list_arg:
                if arg.upper() in map(str.upper, dict_parameters.keys()):
                        string += "${arg.upper()} "
                else:
                        string += arg+" "

        return(string)

                


SCRIPTDIR = os.path.join(os.getcwd(),'resources/submission')

dict_parameters = {}

is_poly = True

file = open(os.path.join(XnfsDir, os.path.join(outputDir,"input_slurm.cfg")),'r')
for line in file.readlines():
        if line.split(' = ')[0] == "Nstat":
                dict_parameters["N"] = [str(i) for i in range(int(line.split(' = ')[1]))]
        else:
                dict_parameters[line.split(' = ')[0]] = string_to_list(line.split(' = ')[1])
                if line.split(' = ')[1].strip() == 'data/toy_domain.in':
                        is_poly = False
file.close()


nb_job = int(np.prod(np.array([len(values) for values in dict_parameters.values()])))

nb_batch = nb_job//96
if nb_job%96 != 0:
        raise IOError("Number of job (%d) is not a multiple of 96" % nb_job)


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


# Max. walltime
WTIME = "1-00:00:00"

# Partition
QUEUE = "Cascade"

# Max. memory per task
MAXMEM = "2G"

# sbatch arguments


for batch in range(nb_batch):


        # Job Name
        JOBNAME = ""
        for keys in dict_mapping_parameters.keys():
                JOBNAME += (keys+"_")
        JOBNAME += str(batch)

        compteur_batch += 1

        file = open(f"resources/submission/tmp/slurm_sweep_process_{outputDir_format}_{batch}_{compteur_batch}.sh", "w")

        file.write("#!/bin/bash\n##\n##  slurm_sweep_tmp.sh\n##  LatticePoly\n")
        file.write(f"##\n##  Created by ppuel on {time.localtime().tm_mday}/{time.localtime().tm_mon}/{time.localtime().tm_year}\n")
        file.write(f"##  Copyright © {time.localtime().tm_year} ENS Lyon. All rights reserved.\n##\n")

        file.write("#SBATCH -o %A_%a.out\n#SBATCH -e %A_%a.err\n#\n")

        file.write(f"#SBATCH --job-name={JOBNAME}        # job name\n")
        file.write(f"#SBATCH --partition={QUEUE}         # partition\n")
        file.write(f"#SBATCH --cpus-per-task=1           # 1 CPU per task\n")
        file.write(f"#SBATCH --mem-per-cpu={MAXMEM}      # 1GiB by CPU\n")
        file.write(f"#SBATCH --ntasks=96                 # 1 tasks\n")
        file.write(f"#SBATCH --time={WTIME}              # one day max\n\n")

        file.write("# Output directory\n")
        file.write(f"OUTPUTDIR={outputDir}\n")

        file.write("\n# Script (relative) path\n")
        file.write(f"SCRIPTDIR={SCRIPTDIR}\n")

        file.write("\n# Data directory\n")
        file.write("XNFSDIR=/Xnfs/physbiochrom/${LOGNAME}/data\n")

        file.write("\n# Error directory\n")
        file.write("ERRORDIR=/Xnfs/physbiochrom/${LOGNAME}/data/${OUTPUTDIR}/tmp\n")

        file.write("# Create Output directory if necessary\n[ ! -d \"${ERRORDIR}\" ] && mkdir -p ${ERRORDIR}\n\n")

        file.write("\n# Scratch directory\n")
        file.write("SCRATCHDIR=/scratch/Cascade/${LOGNAME}/data\n")

        file.write("\n# Relative path to code root directory\n")
        file.write("ROOTDIR=${SCRIPTDIR}/../..\n\n")

        file.write("# Set working directory to root\ncd ${ROOTDIR}\n\n")
        file.write("LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTDIR/lib\n")
        file.write("PYTHONPATH=$PYTHONPATH:$ROOTDIR/resources/h5py\n")
        file.write("export LD_LIBRARY_PATH\nexport PYTHONPATH\n\n")

        file.write("# Values of the parameters\n")

        for keys, values in dict_parameters.items():
                if len(values) == 1:
                        file.write(f"{keys.upper()}={values[0]}\n")
                else:
                        file.write(f"\n{keys.upper()}_ARRAY=({' '.join(dict_parameters_for_job[keys][batch*96:(batch+1)*96])})\n")
                        file.write(f"{keys.upper()}=$"+"{"+f"{keys.upper()}_ARRAY[$"+"[$SLURM_ARRAY_TASK_ID-1]]}\n\n")

        file.write("# Input directory on Xnfs\n")

        input_path = "${XNFSDIR}/${OUTPUTDIR}"
        for keys, values in dict_mapping_parameters.items():
                if keys != "N":
                        input_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
        input_path += "/N/${N}"

        file.write(f"INDIR={input_path}\n\n")

        file.write("# Output directory on scratch\n")

        output_path = "${XNFSDIR}/${OUTPUTDIR}"
        for keys, values in dict_mapping_parameters.items():
                if keys != "N":
                        output_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
        output_path += "/N/${N}"

        file.write(f"OUTDIR={output_path}\n\n")

        file.write("# Create Output directory if necessary\n[ ! -d \"${OUTDIR}\" ] && mkdir -p ${OUTDIR}\n\n")

        if rm_last_process:
                file.write("\n# Remove precedent process data\n")
                file.write("rm -f ${OUTDIR}/process.h5\n")
                
        file.write("\n# Perform processing analysis\n")

        for program, arg_list in dict_program.items():
                file.write(f".venv_bis/bin/python3 {program}"+" ${INDIR}/ ${OUTDIR}/ "+arg_to_variable_slurm(arg_list, dict_parameters)+">> ${OUTDIR}/"+f"process_{program.split('/')[-1].split('.')[0]}_"+"${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out\n")
                
        file.write("\n# Move slurm error/out files\n")        

        file.write("mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${OUTDIR}/process_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out\n")
        file.write("mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${OUTDIR}/process_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err\n")
        
        file.close()

        # sbatch command
        command = f"sbatch -J {JOBNAME} {SCRIPTDIR}/tmp/slurm_sweep_process_{outputDir_format}_{batch}_{compteur_batch}.sh"
        subprocess.run(command, shell=True, executable="/bin/bash")
