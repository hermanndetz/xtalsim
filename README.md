# xtalsim

Open source suite to simulate III-V semiconductor structures on an atomic scale

In a nutshell it internally creates the desired atomic lattice and enables several manipulations such as automatic creation of energetical optimal interface roughnesses or relaxation. In addition several ways to evaluate and properly display the achieved results are provided.

_Features at a glance:_

* creation of atomic lattices with configurable amount of layers, materials
and thicknesses
* periodic boundary conditions in all three dimensions
* determination of energetically best interface configuration
* energy guided optimization
* evaluation of achieved structures
* graphical representation of achieved data
* precise controllability from config files and command line
* various in- and output formats

The code is written in C++ and requires the standard c++11 or newer. To display
data we use [VTK](http://www.vtk.org), which is however not mandatory
for compiling and executing the code.

## Requirements

xtalsim requires the following software packages for compilation. The software is developed with the given version numbers. Other versions and compilers may work, but are not supported. If you spot some troubles, you're encouraged to open an issue or submit a pull request via GitHub.

* GCC C/C++ compiler (>= 7.3.0)
* CMake (>= 3.10.2)

Optional:
* for visualization: VTK (6.3.0)
* for code documentation: Doxygen (>= 1.8.13)

## Installation

The code can be downloaded from GitHub:

```sh
git clone https://github.com/hermanndetz/xtalsim
```

The following commands can be used to compile the executables:

```sh
mkdir build && cd build
cmake ..
make
make docs # optional, requires doxygen
```

## Getting started

## Publications

You can find a list of publications related to xtalsim (or one of its pre-versions) [here](https://github.com/hermanndetz/xtalsim/wiki/publications).

## Acknowledgements


Development of the initial version of xtalsim by Hermann Detz and Juergen Maier was funded by the Austrian Science Fund (FWF): P26100-N27.
