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
import importlib

utils = importlib.import_module("/home/ppuel/Simulation/LatticePoly/LatticePoly/resources/h5py/")
find_exp = utils.find_exp
exp_mapping = utils.exp_mapping

if __name__ == "__main__":
        if len(sys.argv) < 4:
                print("\033[1;31mUsage is %s experience reference script_1 script_2 ...\033[0m" % sys.argv[0])
                sys.exit()
        else:
                experience = sys.argv[1]
                reference = sys.argv[2]
                list_script = sys.argv[3:]
                
        XnfsDir = "/Xnfs/physbiochrom/ppuel/data/"

        for script in list_script:
                if not os.path.isfile(script):
                        print("\033[1;31mThe script %s is not here\033[0m" % script)
                        sys.exit()
        

        expName = find_exp(experience, XnfsDir)

        if reference != "-1":
                refName = find_exp(reference, XnfsDir)
        
        expName_format = expName.replace("/","_")

expDir = os.path.join(XnfsDir, expName)

if reference != "-1":
        refDir = os.path.join("/scratch/Cascade/ppuel/data/", refName)
        refDir_on_Xnfs = os.path.join(XnfsDir, refName)
        

scratchDir = os.path.join("/scratch/Cascade/ppuel/data/", expName)

if input("HP1 ? [Yes]").lower() == 'yes':
        tmp_path = "include/globals.hpp"
else:
        tmp_path = "../LatticePoly_Lucy/include/globals.hpp"


with open(os.path.join(os.getcwd(),tmp_path)) as globalfile:
        for line in globalfile.readlines():
                if "#define L " in line:
                        meta_parameter_L = line.split("#define L ")[1].strip()



scriptDir = os.path.join(os.getcwd(),'resources/submission')

dict_parameters, is_poly, meta_parameter_N, meta_parameter_Nmeas = exp_mapping(expDir)

if reference != "-1":
        ref_parameters, _, ref_parameters_N, ref_parameter_Nmeas = exp_mapping(refDir_on_Xnfs)

        if ref_parameters_Nmeas != meta_parameter_Nmeas:
                print(f"Warning, ref and exp don't have the same Nmeas ({meta_parameter_Nmeas}, {ref_parameters_Nmeas})")
                if input('Yes ? : ').lower()!="yes":
                        sys.exit()

nb_job = int(np.prod(np.array([len(values) for values in dict_parameters.values()])))

job_per_batch = 1

if nb_job%3 == 0:
        job_per_batch = 3

n = nb_job

while n%2==0:
        job_per_batch *= 2
        n //= 2
        if job_per_batch in [32,48]:
                n = 1

nb_batch = nb_job // job_per_batch

print(f"{nb_batch} batchs of {job_per_batch} jobs\n")
if input('Yes ? : ').lower()!="yes":
        sys.exit()

dict_mapping_parameters = {}

for keys,values in dict_parameters.items():
        if len(values) != 1:
                dict_mapping_parameters[keys] = len(values)

if reference != "-1":        
        ref_mapping_parameters = {}

        for keys,values in ref_parameters.items():
                if len(values) != 1:
                        ref_mapping_parameters[keys] = len(values)


dict_parameters_for_job = {keys : [] for keys in dict_mapping_parameters.keys()}

if reference != "-1":
        ref_parameters_for_job = {}
        for keys,values in ref_mapping_parameters.items():
                if values != 1 and keys in dict_parameters.keys():
                        ref_parameters_for_job[keys] = [] 
                        if not set(dict_parameters[keys]).issubset(ref_parameters[keys]):
                                print(f"values of {keys} are not find on ref\nexp : {dict_parameters[keys]}\nref : {ref_parameters[keys]}")
                                sys.exit()


compteur = 1

for keys, values in dict_mapping_parameters.items():
        for job in range(nb_job):
                dict_parameters_for_job[keys].append(dict_parameters[keys][(job//compteur)%values])
                if reference != "-1" and keys in ref_parameters.keys() and len(ref_parameters[keys]) != 1:
                        ref_parameters_for_job[keys].append(dict_parameters[keys][(job//compteur)%values])

        compteur *= values



if reference != "-1":
        for keys, values in ref_parameters_for_job.items():
                if len(values) == 0:
                        ref_parameters_for_job[keys] = [dict_parameters[keys][0] for _ in range(nb_job)]

        print(f"reference parameter {ref_parameters_for_job}")
        if input('Yes ? : ').lower()!="yes":
                sys.exit()


compteur_batch = 0

for batch in range(nb_batch):

        compteur_batch += 1

        file = open(f"resources/submission/tmp/slurm_sweep_process_{expName_format}_{batch}_{compteur_batch}.sh", "w")

        file.write("#!/bin/bash\n##\n##  slurm_sweep_tmp.sh\n##  LatticePoly\n")
        file.write(f"##\n##  Created by ppuel on {time.localtime().tm_mday}/{time.localtime().tm_mon}/{time.localtime().tm_year}\n")
        file.write(f"##  Copyright © {time.localtime().tm_year} ENS Lyon. All rights reserved.\n##\n\n")

        file.write("#SBATCH -o %A_%a.out\n#SBATCH -e %A_%a.err\n\n")

        file.write("\n# Relative path to code root directory\n")
        file.write(f"ROOTDIR={os.getcwd()}\n")

        file.write("\n# Input directory\n")
        file.write(f"EXPDIR={expDir}\n")

        file.write("\n# Output directory\n")
        file.write(f"SCRATCHDIR={scratchDir}\n\n")

        file.write("# Set working directory to root\ncd ${ROOTDIR}\n\n")
        file.write("LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTDIR/lib\n")
        file.write("PYTHONPATH=$PYTHONPATH:$ROOTDIR/resources/h5py\n")
        file.write("export LD_LIBRARY_PATH\nexport PYTHONPATH\n\n")

        file.write("# Values of the parameters\n")

        for keys, values in dict_parameters.items():
                if len(values) == 1:
                        file.write(f"{keys.upper()}={values[0]}\n")
                else:
                        file.write(f"\n{keys.upper()}_ARRAY=({' '.join(dict_parameters_for_job[keys][batch*job_per_batch:(batch+1)*job_per_batch])})\n")
                        file.write(keys.upper()+"=${"+keys.upper()+"_ARRAY[$[$SLURM_ARRAY_TASK_ID-1]]}\n\n")
        
        if reference != "-1":
                for keys, values in ref_parameters.items():
                        if len(values) != 1:
                                file.write(f"\n{keys.upper()}_REF_ARRAY=({' '.join(ref_parameters_for_job[keys][batch*job_per_batch:(batch+1)*job_per_batch])})\n")
                                file.write(keys.upper()+"_REF=${"+keys.upper()+"_REF_ARRAY[$[$SLURM_ARRAY_TASK_ID-1]]}\n\n")

        file.write("\n# Input directory on Xnfs\n")

        input_path = "${EXPDIR}"
        for keys, values in dict_mapping_parameters.items():
                input_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
        
        file.write(f"INPUTDIR={input_path}\n\n")

        if reference != "-1":
                file.write("\n# Reference directory on scratch\n")

                ref_path = refDir
                for keys, values in ref_mapping_parameters.items():
                        ref_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"
                
                file.write(f"REFDIR={ref_path}\n\n")
        else: 
                file.write("\n# No Reference\n")

                file.write(f"REFDIR=-1\n\n")


        file.write("# Output directory on scratch\n")

        output_path = "${SCRATCHDIR}"
        for keys, values in dict_mapping_parameters.items():
                output_path += f"/{keys.upper()}/$"+"{"+f"{keys.upper()}"+"}"

        file.write(f"OUTPUTDIR={output_path}\n\n")

        file.write("# Create output directory if necessary\n[ ! -d \"${OUTPUTDIR}\" ] && mkdir -p ${OUTPUTDIR}\n\n")

        file.write("# Perform post-processing analyses\n")
        for script in list_script:
                file.write(".venv/bin/python3 "+f"{script}"+" ${INPUTDIR}/ ${REFDIR}/ ${OUTPUTDIR}/ "+f"{meta_parameter_N} "+f"{meta_parameter_Nmeas} "+f"{meta_parameter_L} "+"${LDENS} "+f"{is_poly} "+">> ${OUTPUTDIR}/process.out\n")
        
        file.write("\nmv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${OUTPUTDIR}/process_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_.out\n")
        file.write("mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${OUTPUTDIR}/process_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_.err\n")
        

        file.close()

        # Max. walltime
        WTIME = "8-00:00:00"

        # Partition
        QUEUE = "Cascade"

        # Max. memory per task
        MAXMEM = "1G"

        # sbatch arguments
        QARGS=f"--ntasks=1 --mem={MAXMEM} -a 1-{job_per_batch} -t {WTIME} -p {QUEUE}"

        # Job Name
        JOBNAME = ""
        for keys in dict_mapping_parameters.keys():
                JOBNAME += (keys+"_")
        JOBNAME += str(batch)

        # sbatch command
        command = f"sbatch {QARGS} -J {JOBNAME} {scriptDir}/tmp/slurm_sweep_process_{expName_format}_{batch}_{compteur_batch}.sh"
        subprocess.run(command, shell=True, executable="/bin/bash")
        
