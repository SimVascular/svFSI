## Introduction

`svFSI` is a  multi-physics finite element solver designed for computational modeling of the cardiovascular system. Some of the unique capabilities of `svFSI` include modeling cardiac electro-physiology, biological tissue mechanics, blood flow, and large deformation fluid-structure interaction (FSI). `svFSI` also offers a wide choice of boundary conditions for performing patient-specific modeling of cardiovascular biomechanics. The code is parallelized using message-passing-interface (MPI) and offers multiple options to choose a linear solver and preconditioner. `svFSI` can be used as part of the [SimVascular](https://simvascular.github.io) software or can be used a stand-alone solver. The solver is released under the MIT License open source initiative.

## Dependence

The following packages are required to build and use `svFSI`.
   - cmake
   - cmake-curses-gui
   - cmake-gui
   - gcc (version>=4.8.5) with gfortran
   - openmpi or mpich
   - blas & lapack
   - trilinos (optional)

On Ubuntu, most of the dependencies can be downloaded using `apt-get`. On Mac OSX, you may use `brew` to install the dependencies.

## Quick Build

Precompiled binaries for Ubuntu and MacOS are available for download from [SimTK](https://simtk.org/frs/index.php?group_id=188). 

Users are recommended to build from the source code to access the most recent features and bug fixes. A short build instruction is provided here for Linux system.

1. Clone or download the current repository.
2. Create a `build` directory
   ```bash
   cd svFSI && mkdir build && cd build 
   ```
3. Initiate the CMake terminal interface to generate makefiles.
   ```bash
   ccmake ..
   ```
4. This will automatically search for compilers. Follow instructions if necessary. Press “c” to configure repeatedly until CMake parameters no longer change and CMake presents the option “g” for generation. Press “g” to create makefiles and exit. Run `make` in the build directory:
   ```bash
   make 
   ```
   Successful build will generate a solver binary, called `svFSI` in the following directory `build/svFSI-build/bin`.

   For more advanced users, please refer INSTALL.md for detailed platform-specific instructions to install `svFSI`.

## Build With Trilinos

`svFSI` also supports build with [Trilinos](https://github.com/trilinos/Trilinos). Users can build Trilinos locally following its [online documentation](https://docs.trilinos.org/files/TrilinosBuildReference.html). 

The recommended Trilinos third-party libraries (TPLs) include Boost, BLAS, HDF5, HYPRE, LAPACK, and MPI. The required Trilinos packages are Amesos, AztecOO, Epetra, EpetraEXT, Ifpack, ML, MueLU, ROL, Sacado, Teuchos, and Zoltan. 

To enable Trilinos in `svFSI`, users need to turn on the option `SV_USE_TRILINOS` located in the file [`Code/CMake/SimVascularOptions.cmake`](./Code/CMake/SimVascularOptions.cmake) as,

```bash
option(SV_USE_TRILINOS "Use Trilinos Library with svFSI" ON)
``` 

In most cases, users can proceed to build `svFSI` following the [Quick Build](#quick-build), and CMake should be able to locate Trilinos automatically through `find_package`. In case the automatic way fails, users can also specify the path to Trilinos through `ccmake -DCMAKE_PREFIX_PATH:PATH="<Path to Trilinos>/lib/cmake/Trilinos;<Path to other package"`.

For more detailed instructions, please refer INSTALL.md.

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
