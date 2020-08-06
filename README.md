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


## Run

To compile and execute the code, set the simulation parameters to their desired values in the `data/input.cfg` file and type:

~~~shell
make run
~~~

Note that the polymer chain and lattice dimensions are currently hard-coded in the `include/globals.hpp` header to exploit the performance gains of static arrays. Changing their respective values thus requires a full code recompilation.


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
