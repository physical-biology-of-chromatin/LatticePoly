##
##  submit_slurm_process.py
##  LatticePoly
##
##  Created by ppuel on 29/10/2024.
##  Copyright © 2024 ENS Lyon. All rights reserved.
##

import os
import subprocess
import sys
import time

import numpy as np

if len(sys.argv) < 3:
    print(
        "\033[1;31mUsage is %s experience SimProgram1 data1 SimProgram2 data2 ...\033[0m"
        % sys.argv[0]
    )
    sys.exit()
else:
    experience = sys.argv[1]

    dict_program: dict[str, list[str]] = {}

    arg = 2
    tmp_arg = -1
    while arg < len(sys.argv):
        if "resources/" in sys.argv[arg]:
            dict_program[sys.argv[arg]] = []
            tmp_arg = arg
            arg += 1
        else:
            dict_program[sys.argv[tmp_arg]].append(sys.argv[arg])
            arg += 1

    XnfsDir = "/Xnfs/physbiochrom/ppuel/data/"
    inputDir = ""

    verif = 0

    for tmp_experience in os.listdir(XnfsDir):
        if tmp_experience.split("_")[0] == f"EXP{experience}":
            inputDir = tmp_experience
            verif += 1

    if verif != 1:
        print(f"The experience {experience} is not found or found multiple time")
        sys.exit()

    inpurDir_format = inputDir.replace("/", "_")


def string_to_list(input):
    if "," not in input:
        return [input.strip()]
    else:
        return [i.strip() for i in input.split(",")]


def arg_to_variable_slurm(dict_program, dict_parameters):
    string_param_to_python = ""
    list_param_to_python = []
    dict_param_to_python = {}

    tmp_c = 2
    for progId, program in enumerate(dict_program.keys()):
        dict_param_to_python[program] = "+ f'" if len(dict_program[program]) > 0 else ""
        for argId, arg in enumerate(dict_program[program]):
            if arg.upper() in map(str.upper, dict_parameters.keys()):
                string_param_to_python += f"${arg.upper()} "
            else:
                string_param_to_python += f"{arg} "

            list_param_to_python.append(f"P{progId}_{argId} = sys.argv[{tmp_c}]\n")
            dict_param_to_python[program] += "{" + f"P{progId}_{argId}" + "} "
            tmp_c += 1
        dict_param_to_python[program] += "'" if len(dict_program[program]) > 0 else ""

    return (string_param_to_python, list_param_to_python, dict_param_to_python)


SCRIPTDIR = os.path.join(os.getcwd(), "resources/submission")

dict_parameters = {}

is_poly = True

ordinate_list_parameter = []
ordinate_list_length_values = []

file = open(os.path.join(XnfsDir, os.path.join(inputDir, "input_slurm.cfg")), "r")
for line in file.readlines():
    if line.split(" = ")[0] == "Nstat":
        metaParameterN = int(line.split(" = ")[1])
    elif line.split(" = ")[1][0] == "*":
        pass
    else:
        tmp_list_values = string_to_list(line.split(" = ")[1])
        dict_parameters[line.split(" = ")[0]] = tmp_list_values
        if line.split(" = ")[1].strip() == "data/toy_domain.in":
            is_poly = False
        if len(tmp_list_values) > 1:
            ordinate_list_parameter.append(line.split(" = ")[0])
            ordinate_list_length_values.append(len(tmp_list_values))
file.close()


nb_task = int(np.prod(np.array([len(values) for values in dict_parameters.values()])))

MODULO = 1

if nb_task % 3 == 0:
    MODULO = 3

n = nb_task

while n % 2 == 0:
    MODULO *= 2
    n //= 2
    if MODULO in [64, 48]:
        n = 1

nb_job = nb_task // MODULO

NB_NODE_IN_PARALLEL = 16


dict_parameters_for_task = {parameter: [] for parameter in ordinate_list_parameter}

compteur = 1

for paramId, parameter in enumerate(ordinate_list_parameter):
    tmp_length = ordinate_list_length_values[paramId]
    for task in range(nb_task):
        dict_parameters_for_task[parameter].append(
            dict_parameters[parameter][(task // compteur) % tmp_length]
        )
    compteur *= tmp_length

# Max. walltime
WTIME = "6-00:00:00"

# Partition
QUEUE = "Cascade"

# Max. memory per task
MAXMEM = "1G"

# Job Name
JOBNAME = ""
for keys in dict_program.keys():
    JOBNAME += keys.split("/")[-1].split(".")[0] + "_"
JOBNAME += experience

string_param_to_python, list_param_to_python, dict_param_to_python = (
    arg_to_variable_slurm(dict_program, dict_parameters)
)


with open(
    f"resources/submission/tmp/which_python_SimProcess_{JOBNAME}.py", "w"
) as wfile:
    wfile.write(
        "#!/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3\n##\n##  which_python_SimProcess.py\n##  LatticePoly\n"
    )
    wfile.write(
        f"##\n##  Created by ppuel on {time.localtime().tm_mday}/{time.localtime().tm_mon}/{time.localtime().tm_year}\n"
    )
    wfile.write(
        f"##  Copyright © {time.localtime().tm_year} ENS Lyon. All rights reserved.\n##\n\n"
    )
    wfile.write("import sys, subprocess\n\n")
    wfile.write("INDIR = sys.argv[1]\n")

    for param_string in list_param_to_python:
        wfile.write(param_string)

    wfile.write("\n")

    for program, arg_list in dict_program.items():
        wfile.write(
            f"subprocess.run('/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3 {program}"
            + "' + f' {INDIR} ' "
            + dict_param_to_python[program]
            + ", shell = True, executable = '/bin/bash')\n"
        )

subprocess.run(
    f"chmod u+x resources/submission/tmp/which_python_SimProcess_{JOBNAME}.py",
    shell=True,
    executable="/bin/bash",
)

with open(f"resources/submission/tmp/slurm_sweep_SimProcess_{JOBNAME}.sh", "w") as file:
    file.write("#!/bin/bash\n##\n##  slurm_sweep_SimProcess.sh\n##  LatticePoly\n")
    file.write(
        f"##\n##  Created by ppuel on {time.localtime().tm_mday}/{time.localtime().tm_mon}/{time.localtime().tm_year}\n"
    )
    file.write(
        f"##  Copyright © {time.localtime().tm_year} ENS Lyon. All rights reserved.\n##\n"
    )

    file.write("#SBATCH -o tmp/%A_%a.out\n")
    file.write("#SBATCH -e tmp/%A_%a.err\n")
    file.write(f"#SBATCH --job-name={JOBNAME}\n")  # job name
    file.write(f"#SBATCH --partition={QUEUE}\n")  # partition
    file.write(
        f"#SBATCH --array=1-{nb_job}%{NB_NODE_IN_PARALLEL}\n"
    )  # an array of nb_job with max NB_NODE_IN_PARALLEL
    file.write(f"#SBATCH --ntasks-per-node={MODULO}\n")  # MODULO task per job
    file.write(f"#SBATCH --cpus-per-task=1\n")  # 1 CPU per task\n")
    file.write(f"#SBATCH --mem-per-cpu={MAXMEM}\n")  # 2GiB by CPU\n")
    file.write(f"#SBATCH --time={WTIME}\n")  # six day max\n\n")
    file.write(f"#SBATCH --mail-type=END,FAIL\n")
    file.write(f"#SBATCH --mail-user=paul-swann.puel@ens-lyon.fr\n")

    file.write("# Input directory\n")
    file.write(f"INPUTDIR={inputDir}\n")

    file.write("\n# Script (relative) path\n")
    file.write(f"SCRIPTDIR={SCRIPTDIR}\n")

    file.write("\n# Data directory\n")
    file.write("XNFSDIR=/Xnfs/physbiochrom/${LOGNAME}/data\n")

    file.write("\n# Error directory\n")
    file.write("ERRORDIR=/Xnfs/physbiochrom/${LOGNAME}/data/${INPUTDIR}/tmp\n")

    file.write(
        '# Create error directory if necessary\n[ ! -d "${ERRORDIR}" ] && mkdir -p ${ERRORDIR}\n\n'
    )

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
            file.write(
                f"\n{keys.upper()}_ARRAY=({' '.join(dict_parameters_for_task[keys])})\n"
            )

    file.write("\n# Begining of the loop\n")

    file.write(f"for ((i = 0 ; i < {MODULO} ; i++)); do\n")

    for keys in ordinate_list_parameter:
        file.write(
            f"   {keys.upper()}=$"
            + "{"
            + f"{keys.upper()}_ARRAY[$"
            + "[($SLURM_ARRAY_TASK_ID-1)*"
            + f"{MODULO}"
            + "+$i]]}\n\n"
        )

    file.write("\t# Input directory on Xnfs\n")

    input_path = "${XNFSDIR}/${INPUTDIR}"
    for keys in ordinate_list_parameter:
        if keys != "N":
            input_path += f"/{keys.upper()}/$" + "{" + f"{keys.upper()}" + "}"

    file.write(f"\tINDIR={input_path}\n\n")

    file.write(
        '\t# Create InputDir/out directory if necessary\n\t[ ! -d "${INPUTDIR}/out" ] && mkdir -p ${INPUTDIR}/out\n\n'
    )

    file.write("\n\t# Perform SimProcessing analysis\n")

    file.write(
        f"\tsrun -n1 resources/submission/tmp/which_python_SimProcess_{JOBNAME}.py"
        + " ${INDIR}/ "
        + string_param_to_python
        + ">> ${INDIR}/out/"
        + f"SimProcess_{JOBNAME}_"
        + "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out &\n"
    )

    file.write("\ndone\n")

    file.write("\nwait\n")

    file.write("\n# Move slurm error/out files\n")

    file.write(
        "mv ${SLURM_SUBMIT_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${ERRORDIR}/SimProcess_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out\n"
    )
    file.write(
        "mv ${SLURM_SUBMIT_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${ERRORDIR}/SimProcess_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err\n"
    )
    file.write("\n# Move slurm error/out files\n")

    file.write(
        "mv ${SLURM_SUBMIT_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${ERRORDIR}/SimProcess_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out\n"
    )
    file.write(
        "mv ${SLURM_SUBMIT_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${ERRORDIR}/SimProcess_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err\n"
    )

# sbatch command
command = f"sbatch -J {JOBNAME} {SCRIPTDIR}/tmp/slurm_sweep_SimProcess_{JOBNAME}.sh"
subprocess.run(command, shell=True, executable="/bin/bash")
