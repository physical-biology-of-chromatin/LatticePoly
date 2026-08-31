# LatticePoly

MC simulations of polymer and liquid-liquid phase separation on a face-centered cubic lattice.


## Requirements

* `gcc` >= 4.9 
* `cmake` >= 3.8


## Checkout and compilation

To download the code and compile its external library dependencies, simply copy-paste the following lines into a terminal:

~~~shell
git clone --recursive https://github.com/physical-biology-of-chromatin/LatticePoly.git
cd LatticePoly/LatticePoly
git checkout Painter_Valency_Main
make libhdf5
~~~

You may need to remove by hand the `VTK` folder from the M. Tortora original version.

following which the code may be compiled as usual,

~~~shell
make
~~~


## Run

To execute the code, set the simulation parameters to their desired values in the `data/input.cfg` file and type:

~~~shell
./bin/lat data/input.cfg
~~~

Note that the lattice dimension is currently hard-coded in the `include/globals.hpp` header to exploit the performance gains of static arrays; changing its value thus requires a full code recompilation. Compilation and execution may be both achieved through the single command

~~~shell
make run
~~~


## Output

The output data is provided in the [**HDF5**](https://www.hdfgroup.org/solutions/hdf5/) file format, which may be easily visualised using the script coded by D. Erba by variety of third-party open-source MD visualization software. The ancient way rely on a script : `toVTK.py` that translate the h5 file into vtk files that where readable by `Paraview` (version = 5.8)


## Utilities

A handful of post-processing utilities (e.g. gyration tensor analysis by singular value decomposition, MSD calculations via Fourier transform, ...) may be found in the `resources/hdf5` folder. To run them, it is recommended to create a virtual python environment, through which necessary packages can be download and scripts can be executed i.e.,

~~~shell
<path_to_python> -m venv .venv
source .venv/bin/activate
<path_to_python> <script_name> <arguments>
~~~

where `<path_to_python>`  is the path to the venv python executable. Calling a given script without any arguments will output its required argument list to the terminal. All the necessary python modules (`hdf5`, `matplotlib`, `numba`, `scipy`, `psutil`, `networkx`) are freely available through the  `pip` package manager of the virtual environment.

~~~shell
pip install -r Requirements.txt
~~~

The last version script are all class oriented, receiving as first argument the folder of the trajectory. They shall all be compatible with `submit_slurm_NewpPocess.py`.

In `resources`, there is several folder of previous scripts by former PhD students and postdoc like A. Z. Abdulla in `az_resources`, Maxime Tortora in `mt_resources` and P.S. Puel in `pp_resources`. 

There is also the last version of the submission script that threefold : 

  `submit_slurm_Trajectory.py` launch an hypercube array of trajectories in the `PSMN cluster` with parameters taken from `input_slurm.cfg`. In this config file, every parameter with multiple comma separated value will be mapped on. The output is a `EXP__` folder with a given id. Then the tree structured folders match the mapped parameter in the order of the config file like `EXP1_name_of_the_experiment/FIRST_PARAMETER_MAPPED/value/SECOND_PARAMETER_MAPPED/value/.../N/n/` with n from 0 to the number of replicas defined by Nstat in the config file. The config file is also copied at the root of the `EXP__` folder. Several Utilities scripts are run on the trajectories and the resulting data are stored in each the `N/n` folder. 

  `submit_slurm_NewProcess.py` launch an array of jobs in the PSMN cluster for each trajectories of a given experiment. Given post-processing utilities scripts (and their parameters), it will run the new post-process scripts on every trajectories and save the data in their corresponding process.h5 file. This allows to redo scripts on a already processed experiment.

  `submit_slurm_SimProcess.py` launch an array of jobs in the PSMN cluster for each configuration of parameters of a given experiment. Given the sim-processing utilities scripts (and their parameters), it will run the scripts on group of replica trajectories and save the data in a post_process.h5 located at N/. Those scripts are for example data aggregators. 


## Credits

Implemented and maintained by [Maxime Tortora](mailto:maxime.tortora@ens-lyon.fr), partly based on Fortran code by Daniel Jost for the basic polymer simulation module. Updated by [Paul-Swann Puel](mailto:paul-swann.puel@ens-lyon.fr).
