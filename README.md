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
make libhdf5
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

The output data is provided in the [**HDF5**](https://www.hdfgroup.org/solutions/hdf5/) file format, which may be easily visualised using the script coded by D. Erba by variety of third-party open-source MD visualization software. The ancient way rely on a script : toVTK.py that translate the h5 file into vtk files that where readable by Paraview (version = 5.8)


## Utilities

A handful of post-processing utilities (e.g. gyration tensor analysis by singular value decomposition, MSD calculations via Fourier transform, ...) may be found in the `resources/hdf5` folder. To run them, it is recommended to create a virtual python environment, through which necessary packages can be download and scripts can be executed i.e.,

~~~shell
<path_to_python> -m venv .venv
source .venv/bin/activaye
<path_to_python> <script_name> <arguments>
~~~

where `<path_to_python>`  is the path to the venv python executable. Calling a given script without any arguments will output its required argument list to the terminal. All the necessary python modules (`hdf5`, `matplotlib`, `numba`, `scipy`, `psutil`, `networkx`) are freely available through the  `pip` package manager of the virtual environment.

~~~shell
pip install -r Requirements.txt
~~~


## Credits

Implemented and maintained by [Maxime Tortora](mailto:maxime.tortora@ens-lyon.fr), partly based on Fortran code by Daniel Jost for the basic polymer simulation module. Updated by [Paul-Swann Puel](mailto:paul-swann.puel@ens-lyon.fr).
