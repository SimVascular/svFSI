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

Prebuild executables for Ubuntu and MacOS are available for download from [SimTK](https://simtk.org/frs/index.php?group_id=188). Users are recommended to build from the source code to access the most recent features and bug fixes. A short build instruction is provided here for Linux system.

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

A master template is provided in the current repository, [svFSI.inp](./svFSI.inp). Users are also recommended to go through the input files in the [examples](https://github.com/SimVascular/svFSI-Tests) and modify them for their needs.

MPI run can be initiated through
   ```bash
   mpiexec -np <number of MPI processes>  <Path to Build>/svFSI-build/bin/svFSI <Path to input file>
   ```
## Features

`svFSI` provides the capability to model a variety of physics, such as elastodyanmics, heat transfer, convection, fluid, structure and electrophysiology. It also provides multiple options for each physics to cater to the users' diverse needs. In the following tables, *Abbreviation* refers to the variable name in the source code, and users can use global search tool such as `grep` to locate the implementations in the code; *Full name* refers to the name of the model; *Input keyword* refers to the phrase in the input file that can invoke such model.

1. Available isochoric constitutive models for the structure equation
   | *Abbreviation* | *Full name*                     | *Input keyword*                                         |
   | -------------- | ------------------------------- | ------------------------------------------------------- |
   | stIso\_StVK    | Saint Venant-Kirchhoff          | "stVK", "stVenantKirchhoff"                             |
   | stIso\_mStVK   | modified Saint Venant-Kirchhoff | "m-stVK", "modified-stVK", "modified-stVenantKirchhoff" |
   | stIso\_nHook   | Neo-Hookean model               | "nHK", "nHK91", "neoHookean", "neoHookeanSimo91"        |
   | stIso\_MR      | Mooney-Rivlin model             | "MR", "Mooney-Rivlin"                                   |
   | stIso\_HGO     | Holzapfel-Gasser-Ogden model    | "HGO"                                                   |
   | stIso\_Gucci   | Guccione model                  | "Guccione", "Gucci"                                     |
   | stIso\_HO      | Holzapfel-Ogden model           | "HO", "Holzapfel"                                       |

2. Available volumetric constitutive models for the structure equation
   | *Abbreviation* | *Full name*         | *Input keyword*                          |
   | -------------- | ------------------- | ---------------------------------------- |
   | stVol\_Quad    | Quadratic model     | "quad", "Quad", "quadratic", "Quadratic" |
   | stVol\_ST91    | Simo-Taylor91 model | "ST91", "Simo-Taylor91"                  |
   | stVol\_M94     | Miehe94 model       | "M94", "Miehe94"                         |

3. Available constitutive models for the fluid equation
   | *Abbreviation*  | *Full name*                          | *Input keyword*                  |
   | --------------- | ------------------------------------ | -------------------------------- |
   | viscType\_Const | Constant viscosity (Newtonian model) | "constant", "const", "newtonian" |
   | viscType\_CY    | Carreau-Yasuda non-Newtonian model   | "carreau-yasuda", "cy"           |
   | viscType\_Cass  | Cassons non-Newtonian model          | "cassons", "cass"                |

4. Available cardiac electrophysiology models
   | *Abbreviation* | *Full name*                      | *Input keyword*               |
   | -------------- | -------------------------------- | ----------------------------- |
   | cepModel\_AP   | Aliev-Panfilov model             | "ap", "aliev-panfilov"        |
   | cepModel\_BO   | Bueno-Orovio-Cherry-Fenton model | "bo", "bueno-orovio"          |
   | cepModel\_FN   | Fitzhugh-Nagumo model            | "fn", "fitzhugh-nagumo"       |
   | cepModel\_TTP  | tenTusscher-Panfilov model       | "ttp", "tentusscher-panfilov" |

## Additional Resource
More details can be found here:
- FSI: https://simvascular.github.io/docssvFSI.html
- Cardiac mechanics modeling:  https://simvascular.github.io/docsCardiacMechanicsModeling.html
- Eletrophysiology: https://simvascular.github.io/docsElectrophysiology.html

## License
`svFSI` is published under BSD 3-Clause License. Details are included in each source code file.

## Citation
In preparation.
