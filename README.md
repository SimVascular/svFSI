## Introduction

`svFSI` is a  multi-physics finite element solver designed for computational modeling of whole heart dynamics. It is capable of blood flow modeling, large-deformation fluid-structure interaction and electrophysiology modeling. `svFSI` can be used as part of the [SimVascular](https://simvascular.github.io) software or can be used a stand-alone solver.

## Dependence

The following packages are required to build and use `svFSI`.
   - cmake
   - cmake-curses-gui
   - cmake-gui
   - gcc (version>=4.8.5) with gfortran
   - openmpi or mpich
   - blas & lapack
   - trilinos (optional)

## Quick Build

Prebuild executable for Linux system is available for download from [SimTK](https://simtk.org/frs/index.php?group_id=188). Users are recommended to build from the source code to access the most recent features and bug fixes. A short build instruction is provided here for Linux system.

1. Clone or download the current repository.
2. Create a Build directory
   ```bash
   cd svFSI && mkdir Build && cd Build 
   ```
3. Initiate the cmake terminal interface to generate makefiles.
   ```bash
   ccmake ..
   ```
4. This will automatically search for compilers. Follow instructions if necessary. Press “c” to configure repeatedly until cmake presents the option “g” for generation. Press “g” to create makefiles and exit. Run make in the Build directory:
   ```bash
   make 
   ```
   Successful build will generate a solver binary, called `svFSI` in the following directory `Build/svFSI-build/bin`.

## Build With Trilinos

`svFSI` also supports build with [Trilinos](https://github.com/trilinos/Trilinos). Users can build Trilinos locally following its [online document](https://docs.trilinos.org/files/TrilinosBuildReference.html), and these packages from Trilinos should be installed (`Epetra`, `AztecOO`, `Amesos`, `ML`, `MueLu` and `Ifpack`) to ensure compatibility with `svFSI`. 

To enable Trilinos in `svFSI`, users need to turn on its compiling option located in line 62 of [`./Code/CMake/SimVascularOptions.cmake`](./Code/CMake/SimVascularOptions.cmake). 

In most cases, users can proceed to build `svFSI` following the [Quick Build](#quick-build), and cmake should be able to locate Trilinos automatically through `find_package`. In case the automatic way fails, users can also specify the path to Trilinos through `ccmake -DCMAKE_PREFIX_PATH="<Path to Trilinos>;<Path to other package"`.

## Run Simulation

`svFSI` requires a plain-text input file to specify simulation parameters. The syntax of the input file can be found [here](https://sites.google.com/site/memt63/tools/MUPFES/mupfes-scripting).

A template is provided in the current repository, [svFSI.dat](./svFSI.dat). Users are also recommended to go through the input files in the [examples](https://github.com/SimVascular/svFSI-Tests) and modify them for their needs.

MPI run can be initiated through
   ```bash
   mpiexec -np <number of MPI processes>  <Path to Build>/svFSI-build/bin/svFSI input.dat
   ```
*The source code in this repository currently uses an input solver file that is incompatible with the format of the solver file produce by the SimVascular **SV FSI** tool. The SimTK svFSI installers are compatible with the  SimVascular **SV FSI** tool format but use an older version of svFSI.*

## Features

`svFSI` provides the capability to model a variety of physics, such as heat transfer, convection, fluid, structure and electrophysiology. It also provides options for each physics to cater to the users' diverse needs.

1. Available isochoric constitutive models for the structure equation
   | *Abbreviation* | *Full name*                     |
   | -------------- | ------------------------------- |
   | stIso\_StVK    | Saint Venant-Kirchhoff          |
   | stIso\_mStVK   | modified Saint Venant-Kirchhoff |
   | stIso\_nHook   | Neo-Hookean model               |
   | stIso\_MR      | Mooney-Rivlin model             |
   | stIso\_HGO     | Holzapfel-Gasser-Ogden model    |
   | stIso\_Gucci   | Guccione model                  |
   | stIso\_HO      | Holzapfel-Ogden model           |

2. Available volumetric constitutive models for the structure equation
   | *Abbreviation* | *Full name*         |
   | -------------- | ------------------- |
   | stVol\_Quad    | Quadratic model     |
   | stVol\_ST91    | Simo-Taylor91 model |
   | stVol\_M94     | Miehe94 model       |

3. Available constitutive models for the fluid equation
   | *Abbreviation*  | *Full name*                          |
   | --------------- | ------------------------------------ |
   | viscType\_Const | Constant viscosity (Newtonian model) |
   | viscType\_CY    | Carreau-Yasuda non-Newtonian model   |
   | viscType\_Cass  | Cassons non-Newtonian model          |

4. Available cardiac electrophysiology models
   | *Abbreviation* | *Full name*                      |
   | -------------- | -------------------------------- |
   | cepModel\_AP   | Aliev-Panfilov model             |
   | cepModel\_BO   | Bueno-Orovio-Cherry-Fenton model |
   | cepModel\_FN   | Fitzhugh-Nagumo model            |
   | cepModel\_TTP  | tenTusscher-Panfilov model       |

## Additional Resource
More details can be found here: https://simvascular.github.io/docssvFSI.html

## License
`svFSI` is published under BSD 3-Clause License. Details are included in each source code file.

## Citation
In preparation.
