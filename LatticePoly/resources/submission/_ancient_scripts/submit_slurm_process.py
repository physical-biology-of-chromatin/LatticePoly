##
##  submit_slurm_process.py
##  LatticePoly
##
##  Created by ppuel on 29/10/2024.
##  Copyright © 2024 ENS Lyon. All rights reserved.
##

from math import exp
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
        
def arg_to_variable_slurm(dict_program, dict_parameters):
        string_param_to_python = ""
        list_param_to_python = []
        dict_param_to_python = {}

        tmp_c = 3
        for progId, program in enumerate(dict_program.keys()):
                dict_param_to_python[program] = "+ f'" if len(dict_program[program]) > 0 else ""
                for argId, arg in enumerate(dict_program[program]):
                        if arg.upper() in map(str.upper, dict_parameters.keys()):
                                string_param_to_python += f"${arg.upper()} "
                        else:
                                string_param_to_python += f"{arg} "
                        
                        list_param_to_python.append(f"P{progId}_{argId} = sys.argv[{tmp_c}]\n")
                        dict_param_to_python[program] += "{"+f"P{progId}_{argId}"+"} "
                        tmp_c += 1
                dict_param_to_python[program] += "'" if len(dict_program[program]) > 0 else ""
                        
        return(string_param_to_python, list_param_to_python, dict_param_to_python)

                
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


nb_task = 6*6 #int(np.prod(np.array([len(values) for values in dict_parameters.values()])))//96

MODULO = 1

if nb_task%3 == 0:
        MODULO = 3

n = nb_task

while n%2==0:
        MODULO *= 2
        n //= 2
        if MODULO in [64,96]:
                n = 1

nb_job = nb_task // MODULO

NB_NODE_IN_PARALLEL = 6




dict_mapping_parameters = {}

for keys,values in dict_parameters.items():
        if len(values) != 1:
                dict_mapping_parameters[keys] = len(values)

dict_parameters_for_task = {keys : [] for keys in dict_mapping_parameters.keys()}
dict_parameters_for_task['N'] = []

compteur = 1

for keys, values in dict_mapping_parameters.items():
        for task in range(nb_task):
                dict_parameters_for_task[keys].append(dict_parameters[keys][(task//compteur)%values])
        compteur *= values

compteur_batch = 0


# Max. walltime
WTIME = "6-00:00:00"

# Partition
QUEUE = "Cascade"

# Max. memory per task
MAXMEM = "2G"

# Job Name
JOBNAME = ""
for keys in dict_program.keys():
        JOBNAME += (keys.split("/")[-1].split('.')[0]+"_")
JOBNAME += experience

string_param_to_python, list_param_to_python, dict_param_to_python = arg_to_variable_slurm(dict_program, dict_parameters)


for keys, values in dict_parameters.items():
                if len(values) > 1:
                        file.write(f"\n{keys.upper()}_ARRAY=({' '.join(dict_parameters_for_task[keys])})\n")

file.write(f"\tINDIR={input_path}\n\n")

        file.write("\t# Output directory on scratch\n")

        output_path = "${XNFSDIR}/${OUTPUTDIR}"
        for keys, values in dict_mapping_parameters.items():
                if keys != "N":
                        output_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
        output_path += "/N/${N}"

        file.write(f"\tOUTDIR={output_path}\n\n")

        file.write("\t# Create Output directory if necessary\n\t[ ! -d \"${OUTDIR}/out\" ] && mkdir -p ${OUTDIR}/out\n\n")

        if rm_last_process:
                file.write("\n\t# Remove precedent process data\n")
                file.write("\trm -f ${OUTDIR}/process.h5\n")
                
        
        file.write("\n\techo ${OUTDIR}\n")

        file.write("\ndone\n")      


with open(f"resources/submission/tmp/which_python_process_{JOBNAME}.py", "w") as wfile:
        wfile.write("#!/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3\n##\n##  which_python_process.py\n##  LatticePoly\n")
        wfile.write(f"##\n##  Created by ppuel on {time.localtime().tm_mday}/{time.localtime().tm_mon}/{time.localtime().tm_year}\n")
        wfile.write(f"##  Copyright © {time.localtime().tm_year} ENS Lyon. All rights reserved.\n##\n\n")
        wfile.write("import sys, subprocess\n\n")
        wfile.write("INDIR = sys.argv[1]\n")
        wfile.write("OUTDIR = sys.argv[2]\n")

        for param_string in list_param_to_python:
                wfile.write(param_string)


        wfile.write("\n")

        for program, arg_list in dict_program.items():
                wfile.write(f"subprocess.run('/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3 {program}"+"' + f' {INDIR} {OUTDIR} ' "+dict_param_to_python[program]+", shell = True, executable = '/bin/bash')\n")
                
subprocess.run(f"chmod u+x resources/submission/tmp/which_python_process_{JOBNAME}.py", shell=True, executable='/bin/bash')

with open(f"resources/submission/tmp/slurm_sweep_process_{JOBNAME}.sh", "w") as file:

        file.write("#!/bin/bash\n##\n##  slurm_sweep_process.sh\n##  LatticePoly\n")
        file.write(f"##\n##  Created by ppuel on {time.localtime().tm_mday}/{time.localtime().tm_mon}/{time.localtime().tm_year}\n")
        file.write(f"##  Copyright © {time.localtime().tm_year} ENS Lyon. All rights reserved.\n##\n")

        file.write("#SBATCH -o tmp/%A_%a_%J_%j_%s.out\n")
        file.write("#SBATCH -e tmp/%A_%a_%J_%j_%s.err\n")
        file.write(f"#SBATCH --job-name={JOBNAME}\n")                           # job name
        file.write(f"#SBATCH --partition={QUEUE}\n")                            # partition
        file.write(f"#SBATCH --array=1-{nb_job}%{NB_NODE_IN_PARALLEL}\n")       # an array of nb_job with max NB_NODE_IN_PARALLEL
        file.write(f"#SBATCH --ntasks-per-node={MODULO}\n")                     # MODULO task per job
        file.write(f"#SBATCH --cpus-per-task=1\n")                              # 1 CPU per task\n")
        file.write(f"#SBATCH --mem-per-cpu={MAXMEM}\n")                         # 2GiB by CPU\n")
        file.write(f"#SBATCH --time={WTIME}\n")                                 # six day max\n\n")

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

        
        file.write("\n# Perform processing analysis\n")

        file.write(f"srun -n{MODULO} -o resources/submission/tmp/which_python_process_{JOBNAME}.py "+" ${INDIR}/ ${OUTDIR}/ "+string_param_to_python+">> ${OUTDIR}/out/"+f"process_{JOBNAME}_"+"${SLURM_ARRAY_JOB_ID}.out\n")
        
        for keys, values in dict_parameters.items():
                if len(values) != 1:
                        file.write(f"   {keys.upper()}=$"+"{"+f"{keys.upper()}_ARRAY[$"+"[($SLURM_ARRAY_TASK_ID-1)*"+f"{MODULO}"+"+$i]]}\n\n")

        file.write("\t# Input directory on Xnfs\n")

        input_path = "${XNFSDIR}/${OUTPUTDIR}"
        for keys, values in dict_mapping_parameters.items():
                if keys != "N":
                        input_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
        input_path += "/N/${N}"

        

        file.write("\n# Move slurm error/out files\n")        

        file.write("mv ${SLURM_SUBMIT_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${ERRORDIR}/process_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out\n")
        file.write("mv ${SLURM_SUBMIT_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${ERRORDIR}/process_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err\n")

# sbatch command
command = f"sbatch -J {JOBNAME} {SCRIPTDIR}/tmp/slurm_sweep_process_{JOBNAME}.sh"
subprocess.run(command, shell=True, executable="/bin/bash")
