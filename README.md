# xtalsim

Open source suite to simulate III-V semiconductor structures on an atomic scale

In a nutshell it internally creates the desired atomic lattice and enables several manipulations such as automatic creation of energetical optimal interface roughnesses or relaxation. In addition several ways to evaluate and properly display the achieved results are provided.

__Features at a glance:__

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

Optional:
CMake includes the file _CMakeUser.txt_, if found in the root directory of xtalsim. This can be used to change the default path of configuration files for the logger.
```python
add_definitions(-DEASYLOGGING_CONF_DIR="path_to_dir_containing_files")
```

## Getting started

The [wiki](https://github.com/hermanndetz/xtalsim/wiki) provides example simulations as well as a detailed description of the individual programs and file formats.

## Publications

You can find a list of publications related to xtalsim (or one of its pre-versions) [here](https://github.com/hermanndetz/xtalsim/wiki/publications).

## Acknowledgements

Development of the initial version of xtalsim by Hermann Detz and Juergen Maier at TU Wien was funded by the Austrian Science Fund (FWF): P26100-N27 and through an APART Fellowship of the Austrian Academy of Sciences. Recent modifications were also performed at the Central European Institute of Technology at Brno University of Technology with support by the ESF under the project CZ.02.2.69/0.0/0.0/16_027/0008371.
