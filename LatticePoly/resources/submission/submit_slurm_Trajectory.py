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

XnfsDir = "/Xnfs/physbiochrom/ppuel/data/"
ROOTDIR = os.getcwd()

os.makedirs(os.path.join("/Xnfs/physbiochrom/ppuel/data/",outputDir), exist_ok=True)

print(subprocess.run(f"cp /home/ppuel/Simulation/LatticePoly/LatticePoly/resources/submission/input_slurm.cfg {os.path.join('/Xnfs/physbiochrom/ppuel/data/',outputDir)}", shell=True, executable="/bin/bash"))


is_poly = True


dict_parameters = {}

ordinate_list_parameter = ["N"] 
ordinate_list_length_values = [] 

file = open(os.path.join(XnfsDir, os.path.join(outputDir,"input_slurm.cfg")),'r')
for line in file.readlines():
        if line.split(' = ')[0] == "Nstat":
                dict_parameters["N"] = [str(i) for i in range(int(line.split(' = ')[1]))]
                ordinate_list_length_values.append(int(line.split(' = ')[1]))
        else:
                tmp_list_values = string_to_list(line.split(' = ')[1])
                dict_parameters[line.split(' = ')[0]] = tmp_list_values
                if line.split(' = ')[1].strip() == 'data/toy_domain.in':
                        is_poly = False
                if len(tmp_list_values)>1:
                        ordinate_list_parameter.append(line.split(' = ')[0])
                        ordinate_list_length_values.append(len(tmp_list_values))
file.close()

ordinate_list_parameter.append(ordinate_list_parameter.pop(0))
ordinate_list_length_values.append(ordinate_list_length_values.pop(0))

nb_task = int(np.prod(np.array([len(values) for values in dict_parameters.values()])))

MODULO = 1

if nb_task%3 == 0:
        MODULO = 3

n = nb_task

while n%2==0:
        MODULO *= 2
        n //= 2
        if MODULO in [32,48]:
                n = 1


nb_job = nb_task // MODULO

NB_NODE_IN_PARALLEL = 16

dict_mapping_parameters = {}

for keys,values in dict_parameters.items():
        if len(values) != 1:
                dict_mapping_parameters[keys] = len(values)


dict_parameters_for_task = {parameter : [] for parameter in ordinate_list_parameter}

compteur = 1

for paramId, parameter in enumerate(ordinate_list_parameter):
        tmp_length = ordinate_list_length_values[paramId]
        for task in range(nb_task):
                dict_parameters_for_task[parameter].append(dict_parameters[parameter][(task//compteur)%tmp_length])
        compteur *= tmp_length

compteur_batch = 0


# Max. walltime
WTIME = "6-00:00:00"

# Partition
QUEUE = "Cascade"

# Max. memory per task
MAXMEM = "1G"

# Job Name
JOBNAME = outputDir_format



with open(f"resources/submission/tmp/slurm_sweep_trajectory_{outputDir_format}.sh", "w") as file:

        file.write("#!/bin/bash\n##\n##  slurm_sweep_tmp.sh\n##  LatticePoly\n")
        file.write(f"##\n##  Created by ppuel on {time.localtime().tm_mday}/{time.localtime().tm_mon}/{time.localtime().tm_year}\n")
        file.write(f"##  Copyright © {time.localtime().tm_year} ENS Lyon. All rights reserved.\n##\n\n")
        file.write("#SBATCH -o tmp/%A_%a.out\n")
        file.write("#SBATCH -e tmp/%A_%a.err\n")
        file.write(f"#SBATCH --job-name={JOBNAME}\n")                           # job name
        file.write(f"#SBATCH --partition={QUEUE}\n")                            # partition
        file.write(f"#SBATCH --array=1-{nb_job}%{NB_NODE_IN_PARALLEL}\n")       # an array of nb_job with max NB_NODE_IN_PARALLEL
        file.write(f"#SBATCH --ntasks-per-node={MODULO}\n")                     # MODULO task per job
        file.write(f"#SBATCH --cpus-per-task=1\n")                              # 1 CPU per task\n")
        file.write(f"#SBATCH --mem-per-cpu={MAXMEM}\n")                         # 2GiB by CPU\n")
        file.write(f"#SBATCH --time={WTIME}\n")                                 # six day max\n\n")
        file.write(f"#SBATCH --mail-type=END,FAIL\n")
        file.write(f"#SBATCH --mail-user=paul-swann.puel@ens-lyon.fr\n")                                 

        file.write("# Output directory\n")
        file.write(f"OUTPUTDIR={outputDir}\n")

        file.write("\n# Temporary directory\n")
        file.write("TEMPORARYDIR=/tmp/${LOGNAME}\n")

        # file.write("\n# Associated scratch directory\n")
        # file.write("SCRATCHDIR=/scratch/Cascade/${LOGNAME}/data\n")

        file.write("\n# Data directory\n")
        file.write("XNFSDIR=/Xnfs/physbiochrom/${LOGNAME}/data\n")
        
        file.write("\n# Error directory\n")
        file.write("ERRORDIR=${XNFSDIR}/${OUTPUTDIR}/tmp\n")

        file.write("# Create error directory if necessary\n[ ! -d \"${ERRORDIR}\" ] && mkdir -p ${ERRORDIR}\n\n")

        file.write("\n# Relative path to code root directory\n")
        file.write(f"ROOTDIR={ROOTDIR}\n\n")

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
                        file.write(f"\n{keys.upper()}_ARRAY=({' '.join(dict_parameters_for_task[keys])})\n")

# First loop


        file.write("\n# Begining of the loop\n")

        file.write(f"for ((i = 0 ; i < {MODULO} ; i++)); do\n") 

        for keys in ordinate_list_parameter:
                file.write(f"\n\t{keys.upper()}=$"+"{"+f"{keys.upper()}_ARRAY[$"+"[($SLURM_ARRAY_TASK_ID-1)*"+f"{MODULO}"+"+$i]]}\n\n")

        file.write("\t# Temporary directory on tmp\n")

        temporary_path = "${TEMPORARYDIR}/${OUTPUTDIR}"
        for keys in ordinate_list_parameter:
                if keys != 'N':
                        temporary_path += f"_{keys.upper()}_$"+"{"+f"{keys.upper()}"+"}"
        temporary_path += "_N_${N}"

        file.write(f"\tTMPDIR={temporary_path}\n\n")

        file.write("\t# Create temporary directory if necessary\n\t[ ! -d \"${TMPDIR}\" ] && mkdir -p ${TMPDIR}\n\n")

        file.write("\n\t# H5 File path\n")
        file.write("\tFILEPATH=${TMPDIR}/traj.h5\n")

        file.write("\t# Substitution strings\n")

        file.write("\tDIRSUB=\"s|\\(outputDir[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${TMPDIR} ;|;\"\n")
        file.write("\tFILSUB=\"s|\\(H5filePath[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${FILEPATH} ;|;\"\n")

        for keys in dict_parameters.keys():
                if not(keys in ["mode", "exponant"]):
                        file.write(f"\t{keys.upper()}SUB=\"s|\\("+f"{keys}"+"[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${"+f"{keys.upper()}"+"} ;|;\"\n")


        file.write("\n\t# Copy input configuration file to output directory, substituting paths and parameter values\n")

        sed_string = "\tsed -e \"${DIRSUB}\"\"${FILSUB}\""
        for keys in dict_parameters.keys():
                if not(keys in ["mode", "exponant"]):
                        sed_string += "\"${"+f"{keys.upper()}SUB"+"}\""

        file.write(sed_string+" < data/input.cfg > ${TMPDIR}/input.cfg\n")

        file.write("\n\t# Run\n\t./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out & id_array[$i]=$!\n\n")



        file.write("\ndone\n")


# Second loop


        file.write(f"for ((i = 0 ; i < {MODULO} ; i++)); do\n\n")

        file.write("\twait ${id_array[$i]}\n\n")

        for keys in ordinate_list_parameter:
                file.write(f"\n\t{keys.upper()}=$"+"{"+f"{keys.upper()}_ARRAY[$"+"[($SLURM_ARRAY_TASK_ID-1)*"+f"{MODULO}"+"+$i]]}\n\n")

         
        file.write("\t# Output directory on Xnfs\n")

        output_path = "${XNFSDIR}/${OUTPUTDIR}"
        for keys in ordinate_list_parameter:
                if keys != "N":
                        output_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
        output_path += "/N/${N}"

        file.write(f"\tOUTDIR={output_path}\n\n")

        file.write("\t# Create Output directory if necessary\n\t[ ! -d \"${OUTDIR}\" ] && mkdir -p ${OUTDIR}\n\n")

        file.write("\t# Temporary directory on tmp\n")

        temporary_path = "${TEMPORARYDIR}/${OUTPUTDIR}"
        for keys in ordinate_list_parameter:
                if keys != 'N':
                        temporary_path += f"_{keys.upper()}_$"+"{"+f"{keys.upper()}"+"}"
        temporary_path += "_N_${N}"

        file.write(f"\tTMPDIR={temporary_path}\n\n")

        file.write("\tcp ${TMPDIR}/traj.h5 ${OUTDIR}/\n")
        file.write("\tcp ${TMPDIR}/input.cfg ${OUTDIR}/\n")
        file.write("\tcp ${TMPDIR}/log.out ${OUTDIR}/\n")
        

        file.write("\n\t# Perform post-processing analyses\n")
        file.write("\t/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3 resources/h5py/Liq_Density.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
        file.write("\t/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3 resources/h5py/Liq_Cluster.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
        file.write("\t/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3 resources/h5py/Liq_MSD.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
        file.write("\t/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3 resources/h5py/Liq_Droplet.py ${TMPDIR} >> ${TMPDIR}/process.out\n")

        if is_poly:
                file.write("\t/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3 resources/h5py/Poly_MSD.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
                file.write("\t/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3 resources/h5py/Poly_Gyration.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
                file.write("\t/home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3 resources/h5py/Liq_Poly_CoM.py ${TMPDIR} >> ${TMPDIR}/process.out\n")
                 
         
        file.write("\n\t# Move slurm error/out files\n")        


        file.write("\n\t# Move all files to XNFS directory\n")
        file.write("\tmv ${TMPDIR}/process.out ${OUTDIR}/\n")
        file.write("\tmv ${TMPDIR}/process.h5 ${OUTDIR}/\n")
        # file.write("\tmv ${TMPDIR}/liq_graph ${OUTDIR}/\n")
        file.write("\tmv ${TMPDIR}/liq_droplets.pickle ${OUTDIR}/\n")
        file.write("\tmv ${TMPDIR}/liq_simple_droplets.pickle ${OUTDIR}/\n")

        file.write("\n\t# Clean scratch\n")
        file.write("\trm -rf ${TMPDIR}\n")


        file.write("\ndone\n")      

        file.write("\nwait\n")

        file.write("mv ${SLURM_SUBMIT_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${ERRORDIR}/process_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out\n")
        file.write("mv ${SLURM_SUBMIT_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${ERRORDIR}/process_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err\n")


# sbatch command
command = f"sbatch -J {JOBNAME} {ROOTDIR}/resources/submission/tmp/slurm_sweep_trajectory_{outputDir_format}.sh"
subprocess.run(command, shell=True, executable="/bin/bash")
