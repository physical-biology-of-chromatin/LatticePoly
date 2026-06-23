##
##  submit_slurm_Trajectory.py
##  LatticePoly
##
##  Created by ppuel on 29/10/2024.
##  Copyright © 2024 ENS Lyon. All rights reserved.
##

import os
import sys
import time
import subprocess

import numpy as np

class SubmitSlurmTrajectory():
    """SubmitSlurmTrajectory is a class for submitting an array of
    trajectories, encapsulated into a folder named
    EXP{number: int}_{comments: str} into the PSMN, using slurm.

        It creates a bash script named slurm_sweep_trajectory_{exp}.sh and
    uses the command sbatch to launch it into the nodes of the cluster.
    Several python process are computed on each trajectories.

        All the trajectories and process are saved into the Xnfs, in a
    data/exp/ folder, with the input for this experiment and the error
    files.

    The comments are organise so that the majority is at the begining
    of the script, to not pollute the code too much. Because __init__()
    as so many attributes, they are subdivided by which are initialized
    by which init_function (those with a "_" at the begining)
    There are also general comments for those subdivisions with an emphasis
    on the "why ?"
    Then the comments of all the methods of the class are presented with
    an emphasis on the "how ?"
    Then there is most of the code. It's important to remember that this
    script writes also the comments of the bash scripts.

    The main goal of this script is to provide a user agnostic tool to
    start simulation in a nice way, without having to fight again with
    Loïs (PSMN guy (actually very nice)), to understand the subtle
    mystic of this cluster. It is all wrapped up in one sbatch call (see
    after), that launch only one array of jobs which is easy to track.
    There can be RAM allocation that is a bit more precise than before,
    and the scratch is not necessary for those simulations (truth from
    Loïs himself). No more crash, unexpected missing file, IO flood,
    scratch unexpected format, and other weird stuffs. Plus you receive
    an nice e-mail when the simulations are finished...
   	What more could you ask for?


    Parameters
    ----------
    exp_name : str
        Must be "EXP{number: int}_{comments: str}" without / if possible.


    General attributes
    __________________
    exp_name : str
        The experiment name without /

    exp_number : str
        The experiment number to order them chronologicaly in
        the data/ folder

    root_dir : str
        The current working directory. It is used to find the user's name.

    log_name : str
        The user's name on the PSMN

    submission_dir : str
        The resources/submission directory inside LatticePoly

    exp_dir : str
        The experiment wrapping folder in the Xnfs

    venv_location : str
        In case we don't share the same virtual python environment
        name.

    bash_script_path : str
        the path of the bash script, which will to send as a batch of jobs

    Input attributes
    ________________
    self.is_poly : bool
        The way the C++ code is structured doesn't allow a simulation
        not to have a polymer. So if only particles are needed, a toy
        polymer is added, not interacting with the particles.

    self.input_parameters : dict[str, list[str]]
        The all dictonnary of parameter from the input_slurm.cfg, indexed by
        the parameter name. Each values is the list of the parameter
        values mapped for this experiment. If the length of the list is one,
        the parameter is not mapped.

    self.ordered_parameter_names : list[str]
        List of the parameters that will be mapped on. The order is
        important to be maintained for the Xnfs experiment file structure.

    self.ordered_parameter_length : list[int]
        List of the length of the parameters that will be mapped on.


    Sbatch attributes
    _________________
    The way slurm is built is that an array of jobs can be send
    to the cluster nodes. This array of jobs is called a batch.
    Each job can have multiple tasks, to which multiple CPU and
    RAM can be allocated.

    For LatticePoly in PSMN, each task will be the simulation of
    a trajectory and it's processing. One CPU will be allocated
    to each task. The Cascade partition used is composed of several
    nodes, with 96 CPU and 4*96GB of RAM. This is defined by class
    constant PARTITION. The amount of RAM per task is define by
    the class variable MAX_MEM and is usually set to 2GB. More
    in-depth computation could be made at this level to automaticaly
    compute MAX_MEM. Allocating more RAM can decrease the priority
    of the batch in the PSMN queue. It is also the case for maximum
    number of node used at the same time, defined by
    N_NODE_USED_IN_PARALLEL. The maximum run time for a job is
    enforced by the PSMN, and is defined by MAX_RUN_TIME. All those
    information are available at
    https://www.ens-lyon.fr/PSMN/Documentation/clusters_usage/.
    The number of task per node is limited to 32 or 48 to avoid
    bad priority is the queue.

    When the batch is sent, the command squeue --me can be used
    to see its progress. The command scancel --me can be used to
    cancel all your jobs. The slurm documentation is available
    at https://slurm.schedmd.com/.

    self.n_tasks : int
        Total number of trajectories to simulate.

    self.n_task_per_node : int
        Number of CPU used per computing node. It is limited to
        32 or 48 to avoid bad priority in the queue.

    self.n_job : int
        Number of node used to compute all the trajectories.


    Parameter allocation attributes
    _______________________________
    It is a bit ugly, but the way found to pass to each task of
    each job its parameters is to create a massive list of each
    parameter in the slurm file. The structure is :

    p1 = 1, 1, 1,  1, 1, 1,   2, 2, 2,  2, 2, 2
    p2 = 1, 1, 1,  2, 2, 2,   1, 1, 1,  2, 2, 2
    p3 = 1, 2, 3,  1, 2, 3,   1, 2, 3,  1, 2, 3

    for p1 in [1, 2], p2 in [1, 2] and p3 in [1, 2, 3]

    self.input_parameters_for_task : dict[str, list[str]]
        This is the dictonnary which return for each mapping
        simulation parameter its list of values, as described
        in the example above.


    Methods
    -------
    _input_slurm()
        Function that copies the input_slurm config file in the Xnfs
        experiment folder. This file is then read to create the
        input attributes of the class. Nstat, the number of trajectories
        made with the same parameters, is a special parameter. Indeed,
        it is not used to modify input.cfg file of the C++ code. Plus,
        it as to be put at the end of ordered_parameter_names list.

    _sbatch_config()
        Function that compute the parameters required to send the batch
        to the PSMN. It mainly consists on computing the valuation 2-adic
        of the number of tasks.

    _parameter_allocation_to_job()
        Function that construct all the list of parameters for all the tasks.
        Why is it this mathematical equation and not another ? Because this
        one works.

    write()
        Function that write the bash script which will be send to the PSNM
        Most of the comments of this function are written in the bash script.
        Plus, the function is split into sections with comments in this script

        So, at the level of the bash script, the slurm arguments are given with
        the #SBATCH --[arg]=[value]. It is necessary that there is no line without
        a # between the begining of the file and the slurm arguments.

        Then there is the two loop structure. It is a way to assure that all
        trajectories start to be computed, and that when one is finished,
        likely the first one that have started, its processing is launched.
        It is done by the "&" at the end of the run command, that make the script
        to continu to be executed until a "wait" command arrives. It may not
        be the most efficient solution, but it works.

        There is also the :
            ($SLURM_ARRAY_TASK_ID-1)*"+f"{self.n_task_per_node}"+"+$i
        It is the formula which compute the index in the list of parameter
        values for execution of the bash script (job (node)) and for each
        iteration of the for loop (task (trajectory)). How so ?
        $SLURM_ARRAY_TASK_ID is the index given by slurm of job in the job
        array. There is no joke. Here is the quote from the slurm webpage:
        ---
            SLURM_ARRAY_TASK_ID
                Job array ID (index) number.
        ---
        So this formula should be read as :
            Job_index * nb_task_per_node + task_index
        To better understand the subtility of the name, go read the full
        sbatch slurm doc. Hope you won't have to.

        The rest of the code comments are in the bash script.

    execute()
        execute the sbatch command with the bash script.


    string_to_list(input=str)
        Static function that transforms the string of parameters
        in the input file into list of parameters.
    """

    def __init__(self, exp_name):
        # maximum number of node used at the same time
        self.N_NODE_USED_IN_PARALLEL = 16

        # Max. RAM memory allocated per task
        self.MAX_MEM = "2G"

        # Max. run time for a job
        self.MAX_RUN_TIME = "6-00:00:00"

        # Partition used
        self.PARTITION = "Cascade"

        self.exp_name = exp_name.replace("/","_")
        self.exp_number = exp_name.split('_')[0].split('EXP')[1]
        self.root_dir = os.getcwd()
        self.log_name = self.root_dir.split("/")[2]
        self.submission_dir = os.path.join(self.root_dir, "resources/submission")
        xnfs_dir = f"/Xnfs/physbiochrom/{self.log_name}/data/"
        self.exp_dir = os.path.join(xnfs_dir, self.exp_name)
        self.venv_location = [dir for dir in os.listdir(self.root_dir) if "env" in dir][0]
        self.bash_script_path = f"{self.submission_dir}/tmp/slurm_sweep_trajectory_{self.exp_name}.sh"

        # Create the folders if necessary
        os.makedirs(os.path.join(self.submission_dir, 'tmp'), exist_ok = True)
        os.makedirs(self.exp_dir, exist_ok = True)

        self.is_poly = True
        self.input_parameters = {}
        self.ordered_parameter_names = []
        self.ordered_parameter_length = []

        self._input_slurm()

        self.n_tasks = 0
        self.n_task_per_node = 1
        self.n_job = 0

        self._sbatch_config()

        self.input_parameters_for_task = {p_name : [] for p_name in self.ordered_parameter_names}

        self._parameter_allocation_to_job()

    def _input_slurm(self):
        subprocess.run(f"cp {self.submission_dir}/input_slurm.cfg {self.exp_dir}", shell=True, executable="/bin/bash")

        with open(os.path.join(self.exp_dir, "input_slurm.cfg"),'r') as cfg_file:
            for line in cfg_file.readlines():
                param_name, param_value = line.split(' = ')

                if param_name == "Nstat":
                    self.input_parameters["N"] = [str(i) for i in range(int(param_value))]
                    n_stat = int(param_value)
                else:
                    tmp_list_values = self.string_to_list(param_value)
                    self.input_parameters[param_name] = tmp_list_values

                    if param_value.strip() == 'data/toy_domain.in':
                        self.is_poly = False

                    if len(tmp_list_values) > 1:
                        self.ordered_parameter_names.append(param_name)
                        self.ordered_parameter_length.append(len(tmp_list_values))

        self.ordered_parameter_names.append('N') # N is not a simulation parameter so it always needs to
        self.ordered_parameter_length.append(n_stat) # be at the end

    def _sbatch_config(self):
        self.n_tasks = int(np.prod(np.asarray(self.ordered_parameter_length)))

        if self.n_tasks % 3 == 0:
            self.n_task_per_node = 3

        tmp_n = self.n_tasks
        while tmp_n % 2 == 0:
            self.n_task_per_node *= 2
            tmp_n //= 2
            if self.n_task_per_node in [32,48]:
                tmp_n = 1

        self.n_job = self.n_tasks // self.n_task_per_node

    def _parameter_allocation_to_job(self):
        C = 1
        for p_index, p_name in enumerate(self.ordered_parameter_names):
            p_length = self.ordered_parameter_length[p_index]
            for task in range(self.n_tasks):
                self.input_parameters_for_task[p_name].append(
                    self.input_parameters[p_name][(task // C) % p_length]
                )
            C *= p_length

    def write(self):

        with open(self.bash_script_path, "w") as file:

            file.write("#!/bin/bash\n##\n##  slurm_sweep_tmp.sh\n##  LatticePoly\n")
            file.write(f"##\n##  Created by {self.log_name} on {time.localtime().tm_mday}/{time.localtime().tm_mon}/{time.localtime().tm_year}\n")
            file.write(f"##  Copyright © {time.localtime().tm_year} ENS Lyon. All rights reserved.\n##\n#\n")

        # sbatch parameters

            file.write("#SBATCH -o tmp/%A_%a.out\n") #localization of output slurm file
            file.write("#SBATCH -e tmp/%A_%a.err\n") #localization of error slurm file
            file.write(f"#SBATCH --job-name={self.exp_name}\n")
            file.write(f"#SBATCH --partition={self.PARTITION}\n")
            file.write(f"#SBATCH --array=1-{self.n_job}%{self.N_NODE_USED_IN_PARALLEL}\n")
            file.write(f"#SBATCH --ntasks-per-node={self.n_task_per_node}\n")
            file.write("#SBATCH --cpus-per-task=1\n")
            file.write(f"#SBATCH --mem-per-cpu={self.MAX_MEM}\n")
            file.write(f"#SBATCH --time={self.MAX_RUN_TIME}\n")

            if "e_mail.txt" in os.listdir(self.submission_dir):
                with open(os.path.join(self.submission_dir, "e_mail.txt"), 'r') as e_mail_file:
                    file.write(f"#SBATCH --mail-type=END,FAIL\n")   # when the all array is finished, either normaly or by an error
                    file.write(f"#SBATCH --mail-user={e_mail_file.read()}\n")    # send an email to :

        # defining the various directory.

            file.write("# Output directory\n")
            file.write(f"EXP_NAME={self.exp_name}\n")

            file.write("\n# Temporary directory\n")
            file.write("TMP_DIR=/tmp/${LOGNAME}/${EXP_NAME}\n")

            file.write("\n# Data directory\n")
            file.write("EXP_DIR=/Xnfs/physbiochrom/${LOGNAME}/data/${EXP_NAME}\n")

            file.write("\n# Error directory\n")
            file.write("ERROR_DIR=${EXP_DIR}/tmp\n")

            file.write("# Create error directory if necessary\n[ ! -d \"${ERROR_DIR}\" ] && mkdir -p ${ERROR_DIR}\n\n")

            file.write("\n# Relative path to code root directory\n")
            file.write(f"ROOT_DIR={self.root_dir}\n\n")

            file.write("\n# Define python executable\n")
            file.write("PYTHON=${ROOT_DIR}/"+f"{self.venv_location}/bin/python3\n\n")
            file.write("# Executable path\nEXEC=bin/lat\n\n")

        # exporting the libraries

            file.write("# Set working directory to root\ncd ${ROOT_DIR}\n\n")
            file.write("LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOT_DIR/lib\n")
            file.write("PYTHONPATH=$PYTHONPATH:$ROOT_DIR/resources/h5py\n")
            file.write("export LD_LIBRARY_PATH\nexport PYTHONPATH\n\n")

        # defining the list of parameter values

            file.write("# Values of the parameters\n")

            for p_name, parameter_value in self.input_parameters.items():
                if len(parameter_value) == 1:
                    file.write(f"{p_name.upper()}={parameter_value[0]}\n")
                else:
                    file.write(f"\n{p_name.upper()}_ARRAY=({' '.join(self.input_parameters_for_task[p_name])})\n")


        # First loop : the trajectories

            file.write("\n# Begining of the loop\n")
            file.write(f"for ((i = 0 ; i < {self.n_task_per_node} ; i++)); do\n")

            # defining the parameters used for the trajectories

            for p_name in self.ordered_parameter_names:
                file.write(f"\n\t{p_name.upper()}=$"+"{"+f"{p_name.upper()}_ARRAY[$"+"[($SLURM_ARRAY_TASK_ID-1)*"+f"{self.n_task_per_node}"+"+$i]]}\n\n")

            # further directories definition

            file.write("\t# Temporary directory on tmp\n")
            temporary_path = "${TMP_DIR}/"
            for p_name in self.ordered_parameter_names:
                temporary_path += f"_{p_name.upper()}_$"+"{"+f"{p_name.upper()}"+"}"
            file.write(f"\tTMP_OUT_DIR={temporary_path}\n\n")
            file.write("\t# Create temporary directory if necessary\n\t[ ! -d \"${TMP_OUT_DIR}\" ] && mkdir -p ${TMP_OUT_DIR}\n")

            file.write("\n\t# H5 File path\n")
            file.write("\tTMP_TRAJ_PATH=${TMP_OUT_DIR}/traj.h5\n\n")

            # defining the string to substitute the parameters into a copy of the input.cfg file

            file.write("\t# Substitution strings\n")
            file.write("\tOUT_DIR_SUB=\"s|\\(outputDir[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${TMP_OUT_DIR} ;|;\"\n")
            file.write("\tTRAJ_PATH_SUB=\"s|\\(H5filePath[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${TMP_TRAJ_PATH} ;|;\"\n")
            for p_name in self.input_parameters.keys():
                file.write(f"\t{p_name.upper()}_SUB=\"s|\\("+f"{p_name}"+"[[:space:]]*=[[:space:]]*\\)\\(.*;\\)|\\1${"+f"{p_name.upper()}"+"} ;|;\"\n")
            file.write("\n\t# Copy input configuration file to output directory, substituting paths and parameter values\n")
            sed_string = "\tsed -e \"${OUT_DIR_SUB}\"\"${TRAJ_PATH_SUB}\""
            for p_name in self.input_parameters.keys():
                sed_string += "\"${"+f"{p_name.upper()}_SUB"+"}\""
            file.write(sed_string+" < data/input.cfg > ${TMP_OUT_DIR}/input.cfg\n")

            # running the simulation

            file.write("\n\t# Run\n\t./${EXEC} ${TMP_OUT_DIR}/input.cfg > ${TMP_OUT_DIR}/log.out & IDARRAY[$i]=$!\n")

            # end of the first loop

            file.write("\ndone\n\n")


        # Second loop : the processing

            file.write(f"for ((i = 0 ; i < {self.n_task_per_node} ; i++)); do\n\n")

            # assure that the trajectories is finished.

            file.write("\t# Wait for the relative trajectory to be computed\n")
            file.write("\twait ${IDARRAY[$i]}\n\n")

            # redefining the parameters

            file.write("\t# Define the input parameters for this trajectory\n")
            for p_name in self.ordered_parameter_names:
                file.write(f"\t{p_name.upper()}=$"+"{"+f"{p_name.upper()}_ARRAY[$"+"[($SLURM_ARRAY_TASK_ID-1)*"+f"{self.n_task_per_node}"+"+$i]]}\n")

            # even more directories definition

            file.write("\n\t# Output directory on Xnfs\n")
            exp_path = "${EXP_DIR}"
            for p_name in self.ordered_parameter_names:
                    exp_path += f"/{p_name.upper()}/$"+"{"+f"{p_name.upper()}"+"}"
            file.write(f"\tEXP_OUT_DIR={exp_path}\n\n")
            file.write("\t# Create Output directory if necessary\n\t[ ! -d \"${EXP_OUT_DIR}\" ] && mkdir -p ${EXP_OUT_DIR}\n\n")
            file.write("\t# Temporary directory on tmp\n")
            temporary_path = "${TMP_DIR}/"
            for p_name in self.ordered_parameter_names:
                if p_name != 'N':
                    temporary_path += f"_{p_name.upper()}_$"+"{"+f"{p_name.upper()}"+"}"
            temporary_path += "_N_${N}"
            file.write(f"\tTMP_OUT_DIR={temporary_path}\n\n")

            # copy of the trajectories is case of an error during processing

            file.write("\t# Copy the trajectory file before doing the processing\n")
            file.write("\t# in case the one of the python file raises an issue\n")
            file.write("\tcp ${TMP_OUT_DIR}/traj.h5 ${EXP_OUT_DIR}/\n")
            file.write("\tcp ${TMP_OUT_DIR}/input.cfg ${EXP_OUT_DIR}/\n")
            file.write("\tcp ${TMP_OUT_DIR}/log.out ${EXP_OUT_DIR}/\n")

            # execute python processing scripts

            file.write("\n\t# Perform processing analyses\n")
            file.write("\t${PYTHON} resources/h5py/Liq_Density.py ${TMP_OUT_DIR} >> ${TMP_OUT_DIR}/process.out\n")
            file.write("\t${PYTHON} resources/h5py/Liq_Cluster.py ${TMP_OUT_DIR} >> ${TMP_OUT_DIR}/process.out\n")
            file.write("\t${PYTHON} resources/h5py/Liq_MSD.py ${TMP_OUT_DIR} >> ${TMP_OUT_DIR}/process.out\n")
            file.write("\t${PYTHON} resources/h5py/Liq_Droplet.py ${TMP_OUT_DIR} >> ${TMP_OUT_DIR}/process.out\n")

            # even more python script if there is a polymer

            if self.is_poly:
                file.write("\t${PYTHON} resources/h5py/Poly_MSD.py ${TMP_OUT_DIR} >> ${TMP_OUT_DIR}/process.out\n")
                file.write("\t${PYTHON} resources/h5py/Poly_Gyration.py ${TMP_OUT_DIR} >> ${TMP_OUT_DIR}/process.out\n")
                file.write("\t${PYTHON} resources/h5py/Liq_Poly_CoM.py ${TMP_OUT_DIR} >> ${TMP_OUT_DIR}/process.out\n")

            # saving the files in Xnfs

            file.write("\n\t# Move all files to XNFS directory\n")
            file.write("\tmv ${TMP_OUT_DIR}/process.out ${EXP_OUT_DIR}/\n")
            file.write("\tmv ${TMP_OUT_DIR}/process.h5 ${EXP_OUT_DIR}/\n")
            file.write("\tmv ${TMP_OUT_DIR}/liq_droplets.pickle ${EXP_OUT_DIR}/\n")
            file.write("\tmv ${TMP_OUT_DIR}/liq_simple_droplets.pickle ${EXP_OUT_DIR}/\n")

            # release the RAM of the node

            file.write("\n\t# Clean the RAM\n")
            file.write("\trm -rf ${TMP_OUT_DIR}\n")

            # end of the second loop

            file.write("\ndone\n")

        # wait that all the trajectories are processed

            file.write("\nwait\n")

        # save the slurm error file (super important) and the output file (always empty)

            file.write("mv ${SLURM_SUBMIT_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${ERROR_DIR}/process_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out\n")
            file.write("mv ${SLURM_SUBMIT_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${ERROR_DIR}/process_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err\n")

        # maybe one day even post-processing script will be added

          # file.write("${PYTHON} resources/submission/submit_slurm_NewProcess.py " + f"{self.exp_number} resources/h5py/SimLiqRadius" )

    def send_the_batch(self):
        command = f"sbatch -J {self.exp_name} {self.bash_script_path}"
        subprocess.run(command, shell=True, executable="/bin/bash")


    @staticmethod
    def string_to_list(input):
        if not ',' in input:
            return([input.strip()])
        else:
            return([i.strip() for i in input.split(",")])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("\033[1;31mUsage is %s exp_name \033[0m" % sys.argv[0])
        sys.exit()

    exp_name = sys.argv[1]

    submit = SubmitSlurmTrajectory(exp_name)
    submit.write()
    submit.send_the_batch()
