## NOTE

This branch is created from master on Jan 4, 2022 and will be referenced in the manuscript to be submitted for the Journal of Open Source Software (JOSS).

## Introduction

`svFSI` is a  multi-physics finite element solver designed for computational modeling of the cardiovascular system. Some of the unique capabilities of `svFSI` include modeling cardiac electrophysiology, biological tissue mechanics, blood flow, and large deformation fluid-structure interaction (FSI). `svFSI` also offers a wide choice of boundary conditions for performing patient-specific modeling of cardiovascular biomechanics. The code is parallelized using message-passing-interface (MPI) and offers multiple options to choose a linear solver and preconditioner. `svFSI` can be used as part of the [SimVascular](https://simvascular.github.io) software or can be used as a stand-alone solver.

## Dependence

The following packages are required to build and use `svFSI`.
   - cmake
   - cmake-curses-gui
   - cmake-gui
   - gcc (version>=4.8.5) with gfortran
   - openmpi or mpich
   - blas & lapack
   - trilinos (optional)

On Ubuntu, most of the dependencies can be installed using `apt install`. On macOS, the dependencies may be installed using `brew`. Apart from GNU compilers, `svFSI` can also be built with Intel oneAPI Toolkits. For more details, please refer to [`INSTALL.md`](./INSTALL.md#Build) and [`INSTALL-DEPS.md`](./INSTALL-DEPS.md#intel-oneapi-toolkitsd).

## Quick Build

Precompiled binaries for Ubuntu and MacOS are available for download from [SimTK](https://simtk.org/frs/index.php?group_id=188).

Users are recommended to build from the source code to access the most recent features and bug fixes. Instructions for a quick build are provided here for a Linux/Mac OS system.

1. Clone or download the current repository.
2. Create a `build` directory
   ```bash
   cd svFSI && mkdir build && cd build
   ```
3. Initiate the CMake terminal interface to generate makefiles.
   ```bash
   ccmake ..
   ```
4. This will automatically search for compilers. Follow instructions if necessary. Press “c” to configure repeatedly until CMake parameters no longer change and CMake presents the option “g” to generate. Press “g” to create makefiles and exit. Run `make` in the build directory:
   ```bash
   make
   ```
   A successful build will generate a solver binary, called `svFSI` in the following directory `build/svFSI-build/bin`.

   For more advanced users, please refer [`INSTALL.md`](./INSTALL.md) for detailed platform-specific instructions to install `svFSI`.

## Build With Trilinos

`svFSI` also supports compilation with [Trilinos](https://github.com/trilinos/Trilinos). Users can build Trilinos locally following its [online documentation](https://docs.trilinos.org/files/TrilinosBuildReference.html).

The recommended Trilinos third-party libraries (TPLs) include Boost, BLAS, HDF5, HYPRE, LAPACK, MPI, and MUMPS. The required Trilinos packages are Amesos, AztecOO, Epetra, EpetraEXT, Ifpack, ML, MueLU, ROL, Sacado, Teuchos, and Zoltan.

To enable Trilinos in `svFSI`, users need to turn on the option `SV_USE_TRILINOS` located in the file [`Code/CMake/SimVascularOptions.cmake`](./Code/CMake/SimVascularOptions.cmake) as,

```bash
option(SV_USE_TRILINOS "Use Trilinos Library with svFSI" ON)
```

In most cases, users can proceed to build `svFSI` following the [Quick Build](#quick-build), and CMake should be able to locate Trilinos automatically through `find_package`. In case the automatic way fails, users can also specify the path to Trilinos through `ccmake -DCMAKE_PREFIX_PATH:PATH="<Path_to_Trilinos>/lib/cmake/Trilinos;<Path_to_any_other_package>;"`.

For more detailed instructions, please refer INSTALL.md.

## Run Simulation

`svFSI` requires a plain-text input file to specify simulation parameters. The syntax of the input file can be found [here](https://sites.google.com/site/memt63/tools/MUPFES/mupfes-scripting).

A master template is provided in the current repository, [svFSI_master.inp](./svFSI_master.inp). Users are also recommended to go through the input files in the [examples](https://github.com/SimVascular/svFSI-Tests) and modify them for their needs.

An MPI-based run can be initiated through
```bash
mpiexec -np <number of MPI processes>  <Path to Build>/svFSI-build/bin/svFSI <Path to input file>
```
## Features

`svFSI` provides the capability to model a variety of physics including unsteady diffusion, linear and nonlinear elastodynamics, convective heat transfer, fluid flows, fluid-structure-interaction (FSI), and cardiac electrophysiology. As the code is modular, the users are provided with a choice to couple these physics depending on their needs. We strongly recommend users to browse through the examples provided in the GitHub repository [svFSI-Tests](https://github.com/SimVascular/svFSI-Tests) to get a detailed insight into the capability of the code.

Below, we provide a list of the available choice of constitutive models for different types of equations being solved. Users are also encouraged to implement new constitutive models. Users may use global search tools such as `grep` to locate the implementations of the available constitutive models in the code using the abbreviated names below.

*Abbreviation* refers to the variable name in the source code; *Full name* refers to the generic name of the model; *Input keyword* refers to the phrase in the input file that can invoke such model.

1. Available isochoric constitutive models for the structure equation
   | *Abbreviation* | *Full name*                        | *Input keyword*                                         |
   | -------------- | ---------------------------------- | ------------------------------------------------------- |
   | stIso\_stVK    | Saint Venant-Kirchhoff             | "stVK", "stVenantKirchhoff"                             |
   | stIso\_mStVK   | modified Saint Venant-Kirchhoff    | "m-stVK", "modified-stVK", "modified-stVenantKirchhoff" |
   | stIso\_nHook   | Neo-Hookean model                  | "nHK", "nHK91", "neoHookean", "neoHookeanSimo91"        |
   | stIso\_MR      | Mooney-Rivlin model                | "MR", "Mooney-Rivlin"                                   |
   | stIso\_HGO_d   | Holzapfel-Gasser-Ogden (decoupled) | "HGO", "HGO-d", HGO-decoupled"                          |
   | stIso\_HGO_ma  | HGO model (modified anisotropy)    | "HGO-ma", "HGO-modified"                                |
   | stIso\_Gucci   | Guccione model                     | "Guccione", "Gucci"                                     |
   | stIso\_HO_d    | Holzapfel-Ogden model (decoupled)  | "HO", "Holzapfel", "HO-decoupled", "HO-d"               |
   | stIso\_HO_ma   | HO model (modified anisotropy)     | "HO-ma", "HO-modified"                                  |

2. Available volumetric constitutive models for the structure equation
   | *Abbreviation* | *Full name*         | *Input keyword*                          |
   | -------------- | ------------------- | ---------------------------------------- |
   | stVol\_Quad    | Quadratic model     | "quad", "Quad", "quadratic", "Quadratic" |
   | stVol\_ST91    | Simo-Taylor91 model | "ST91", "Simo-Taylor91"                  |
   | stVol\_M94     | Miehe94 model       | "M94", "Miehe94"                         |

3. Available constitutive models for the fluid equation
   | *Abbreviation*  | *Full name*                          | *Input keyword*                  |
   | --------------- | ------------------------------------ | -------------------------------- |
   | viscType\_Const | Constant viscosity (Newtonian model) | "Constant", "Const", "Newtonian" |
   | viscType\_CY    | Carreau-Yasuda non-Newtonian model   | "Carreau-Yasuda", "CY"           |
   | viscType\_Cass  | Cassons non-Newtonian model          | "Cassons", "Cass"                |

4. Available cardiac electrophysiology models
   | *Abbreviation* | *Full name*                      | *Input keyword*               |
   | -------------- | -------------------------------- | ----------------------------- |
   | cepModel\_AP   | Aliev-Panfilov model             | "AP", "Aliev-Panfilov"        |
   | cepModel\_BO   | Bueno-Orovio-Cherry-Fenton model | "BO", "Bueno-Orovio"          |
   | cepModel\_FN   | Fitzhugh-Nagumo model            | "FN", "Fitzhugh-Nagumo"       |
   | cepModel\_TTP  | tenTusscher-Panfilov model       | "TTP", "tenTusscher-Panfilov" |

## Additional Resource
More details can be found here:
- Fluid-Structure Interaction (FSI): https://simvascular.github.io/docssvFSI.html
- SimCardio: http://simvascular.github.io/docsSimCardio.html
- Cardiac electrophysiology modeling: http://simvascular.github.io/docsSimCardio.html#cep-modeling
- Cardiac mechanics modeling:  http://simvascular.github.io/docsSimCardio.html#mechanics-modeling
- Prescribed-motion LV modeling: https://simvascular.github.io/docsSimCardio.html#automatic-cardiac-modeling

## Citation
In preparation.
