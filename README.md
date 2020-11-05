# Replication

The following elements were added to the code:

Parameters
-
int Parameter "replstatus" in MCTad: 0 if unreplicated, 1 if fork moving to the right, -1 if fork moving to the left, -2 if the monomer of the original chain has been used as a template for replication, +2 if the monomer is added through the replication process.

int vector "activeforks" in MCPoly: contains all the monomers which are forks 

int "SisterID" in MCTad: since moving (and thus creating a monomer) doesn't follow linear order the parameter saves the template monomer position in the chain. (I'm not completely sure if it could be useful)

void MCReplicPoly::CreateFork()
-
randomly select a position between the two ends of the chain and after a checking the replstatus(=0) of the 3 monomer involved (2 forks and one on top of which I duplicate) initialize a new origin.

void MCReplicPoly::MoveFork(int forkID, int i)
-
After selecting a random element from the vector activeforks at index i, both of the values serve as an input for the function:

If the fork is in the second or penultimate position it cannot move (in order to save the isFork connotation, but this probably needs to change) otherwise starts to move in the direction dictated by the replstatus parameter.

The general procedure:

1)A new monomer is created and pushed in TadConf.

2)Neighbours and bonds parameters of the closest replicated monomer are updated.

3)The properties of the new monomer (regarding his closest replicated neighbour) are updated.

4)the new bond is created and pushed into TopoConf.

5)The Update() function manages the new number of monomers and initializes the new bond at the fork.

6)All parameters regarding the fork dynamics (link and replstatus) are updated.

When two fork meets (two consecutive fork) a similar procedure is used but 2 monomer and one bonds between them are created: to adjust this inconsistency I set a probability of 0.5 to "open" the replication bubble.
When two fork meets (two consecutive fork) a similar procedure is used but 2 monomer and one bond are created: to adjust this inconsistency I set a probability of 0.5 to "open" the replication bubble.

Comments/ Errors:

In the process of moving in the negative direction I wasn't able to create the fork bond first and then update the second one: to solve it I had to always follow the order described and if needed swap the first and second bond/neighourg (l.262) .


void MCSim<lattice, polymer>::Run()
-
At each run a cumulative number of trial moves (forks moves + tad moves) are made, according to the probability to pick one kind or the other (Ntad/Ntad+forks or forks//Ntad+forks).

When selected the fork movement as a probability to be moved according to an input "Replicationrate"
Similarly at each step a fork may be created according to a certain "Originrate" .

Comments/ Errors

The template functions UpdateFork and UpdateCreate used in MCSim work, but probably is not the proper implementation.
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
make libvtk
~~~

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

The output data is provided in the [**VTK**](https://vtk.org) file format, which may be easily visualised using a variety of third-party open-source software (e.g. [**ParaView**](https://www.paraview.org)).


## Utilities

A handful of post-processing utilities (e.g. gyration tensor analysis by singular value decomposition, MSD calculations via Fourier transform, ...) may be found in the `resources` folder. To run them, it is recommended to install the [**miniconda**](https://docs.conda.io/en/latest/miniconda.html) python distribution, through which they may be executed as regular standalone scripts, i.e.,

~~~shell
<path_to_python> <script_name> <arguments>
~~~

where `<path_to_python>`  is the path to the `miniconda` python executable. Calling a given script without any arguments will output its required argument list to the terminal. All the necessary python modules (`vtk`, `matplotlib`, `numba`, `scipy`, `psutil`, `fileseq`) are freely available through the  `pip` package manager, and may be simply installed in the standard fashion,

~~~shell
<path_to_python> -m pip install <module_name>
~~~


## Credits

Implemented and maintained by [Maxime Tortora](mailto:maxime.tortora@ens-lyon.fr), partly based on Fortran code by Daniel Jost for the basic polymer simulation module.
