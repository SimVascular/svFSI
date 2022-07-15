## NOTE

This branch is created from master on Jan 4, 2022 and will be referenced in the manuscript to be submitted for the Journal of Open Source Software (JOSS).

## Introduction

`svFSI` is a multi-physics finite element solver designed for computational modeling of the cardiovascular system. It is a major component of the ongoing SimVascular [**SimCardio**](http://simvascular.github.io/docsSimCardio.html) project that aims to provide the complete pipeline for cardiac modeling, from image segmentation to computational modeling.

Some of the unique capabilities of `svFSI` include modeling cardiac electrophysiology, biological tissue mechanics, blood flow, and large deformation fluid-structure interaction (FSI). `svFSI` also offers a wide choice of boundary conditions for performing patient-specific modeling of cardiovascular biomechanics. The code is parallelized using message-passing-interface (MPI) and offers multiple options to choose a linear solver and preconditioner. `svFSI` can be used as part of the [SimVascular](https://simvascular.github.io) software or can be used as a stand-alone solver. It is distributed under a MIT-like open source license.

## Binary and Container
Precompiled binaries for Ubuntu and MacOS are available for download from [SimTK](https://simtk.org/frs/index.php?group_id=188).

Instructions to build and run `svFSI` in Docker container are provided [here](./Docker/README.md).

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

## Quick Build from Source

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
   A successful build will generate a solver binary called `svFSI` in the following directory `build/svFSI-build/bin`.

   For more advanced users, please refer [`INSTALL.md`](./INSTALL.md) for detailed platform-specific instructions to install `svFSI`.

## Build With Trilinos

`svFSI` also supports compilation with [Trilinos](https://github.com/trilinos/Trilinos). Users can build Trilinos locally following its [online documentation](https://docs.trilinos.org/files/TrilinosBuildReference.html).

The recommended Trilinos third-party libraries (TPLs) include Boost, BLAS, HDF5, HYPRE, LAPACK, MPI, and MUMPS. The required Trilinos packages are Amesos, AztecOO, Epetra, EpetraEXT, Ifpack, ML, MueLU, ROL, Sacado, Teuchos, and Zoltan.

To enable Trilinos in `svFSI`, users need to turn on the option `SV_USE_TRILINOS` located in the file [`Code/CMake/SimVascularOptions.cmake`](./Code/CMake/SimVascularOptions.cmake) as,

```bash
option(SV_USE_TRILINOS "Use Trilinos Library with svFSI" ON)
```

In most cases, users can proceed to build `svFSI` following the [Quick Build](#quick-build), and CMake should be able to locate Trilinos automatically through `find_package`. In case the automatic way fails, users can also specify the path to Trilinos through `ccmake -DCMAKE_PREFIX_PATH:PATH="<Path_to_Trilinos>/lib/cmake/Trilinos;<Path_to_any_other_package>;"`.

For more detailed instructions, please refer to [`INSTALL.md`](./INSTALL.md).

## Run Simulation

`svFSI` requires a plain-text input file to specify simulation parameters. The syntax of the input file can be found [here](http://simvascular.github.io/docssvFSI.html#input).

A master template is provided in the current repository, [svFSI_master.inp](./svFSI_master.inp). Users are also recommended to go through the input files in the [examples](https://github.com/SimVascular/svFSI-Tests) and modify them for their needs.

An MPI-based run can be initiated through
```bash
mpiexec -np <number of MPI processes>  <Path to Build>/svFSI-build/bin/svFSI <Path to input file>
```
## Features

`svFSI` provides the capability to model a variety of physics including unsteady diffusion, linear and nonlinear elastodynamics, convective heat transfer, fluid flows, fluid-structure-interaction (FSI), and cardiac electrophysiology. As the code is modular, the users are provided with a choice to couple these physics depending on their needs. We strongly recommend users to browse through the examples provided in the GitHub repository [svFSI-Tests](https://github.com/SimVascular/svFSI-Tests) to get a detailed insight into the capability of the code. Also, most of the examples contain established simulation results, which users can use to verify the functionality of `svFSI`. Here is a list of the main features of `svFSI`.

   | *Physics Solved*  | *Documentation/Tutorial* | *Examples* |
   |-------------------|--------------------------|------------|
   | Fluid             | [Webpage](http://simvascular.github.io/docssvFSI-Fluid.html)                                                                                                              | [pipe flow with RCR BC](https://github.com/SimVascular/svFSI-Tests/tree/master/04-fluid/01-pipe3D_RCR);<br>[dye transportation](https://github.com/SimVascular/svFSI-Tests/tree/master/04-fluid/02-dye_AD);<br>[GenBC/cplBC](https://github.com/SimVascular/svFSI-Tests/tree/master/04-fluid/04-3D0D-coupling-BC);<br>[Non-Newtonian flow](https://github.com/SimVascular/svFSI-Tests/tree/master/04-fluid/05-nonNewtonian)                                                                                                                                                                                                            |
   | Structure         | [Webpage](http://simvascular.github.io/docssvFSI-Structure.html); [YouTube](https://www.youtube.com/watch?v=Jm3VSi6Aci8&list=PL1CBZ8Wh-xvRnux0eMmbZPbx-C078Qzqu&index=2)   | struct:<br>[block compression](https://github.com/SimVascular/svFSI-Tests/tree/master/05-struct/01-block-compression);<br>[passive inflation of LV model](https://github.com/SimVascular/svFSI-Tests/tree/master/05-struct/02-LV-Guccione-passive)<br>ustruct:<br>[block compression](https://github.com/SimVascular/svFSI-Tests/tree/master/06-ustruct/01-block-compression);<br>[tension of arterial strip](https://github.com/SimVascular/svFSI-Tests/tree/master/06-ustruct/02-tensile-adventitia_HGO);<br>[active inflation of LV model](https://github.com/SimVascular/svFSI-Tests/tree/master/06-ustruct/03-LV-Guccione-active) |
   | Electrophysiology | [Webpage](http://simvascular.github.io/docssvFSI-CEP.html); [YouTube](https://www.youtube.com/watch?v=TCK3SmGwBa8&list=PL1CBZ8Wh-xvRnux0eMmbZPbx-C078Qzqu&index=2) | [Aliev-Panfilov model](https://github.com/SimVascular/svFSI-Tests/tree/master/08-cep/01-2Dsqr_AP); <br>[ten-Tusscher-Panfilov model](https://github.com/SimVascular/svFSI-Tests/tree/master/08-cep/03-benchmark_tTP); <br>[Bueno-Orovio-Cherry-Fenton model](https://github.com/SimVascular/svFSI-Tests/tree/master/08-cep/04-2Dspiral_BO); <br>[Purkinje network](https://github.com/SimVascular/svFSI-Tests/tree/master/08-cep/05-Purkinje)                                                                                                                                                                                          |
   | FSI               | [Webpage](http://simvascular.github.io/docssvFSI-FSI.html); [YouTube](https://www.youtube.com/watch?v=QIpyThIAD7k&list=PL1CBZ8Wh-xvRnux0eMmbZPbx-C078Qzqu&index=4)     | ALE:<br>[2D heart valve](https://github.com/SimVascular/svFSI-Tests/tree/master/07-fsi/ale/01-channel-leaflets_2D); <br>[2D flag behind a block](https://github.com/SimVascular/svFSI-Tests/tree/master/07-fsi/ale/02-channel-block-flag_2D); <br>[pressure pulse inside aorta](https://github.com/SimVascular/svFSI-Tests/tree/master/07-fsi/ale/03-pipe_3D)<br>CMM:<br>[pipe flow with RCR BC](https://github.com/SimVascular/svFSI-Tests/tree/master/07-fsi/cmm/01-pipe_RCR);<br>[vein graft](https://github.com/SimVascular/svFSI-Tests/tree/master/07-fsi/cmm/02-vein-graft)                                                      |


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

## Documentation
More details can be found on the [**svFSI**](http://simvascular.github.io/docssvFSI.html) page, and direct links to the documentation for different functionalities are provided here:
- [Fluid-Structure Interaction (FSI)](http://simvascular.github.io/docssvFSI-FSI.html)
- [Cardiac electrophysiology modeling](http://simvascular.github.io/docssvFSI-CEP.html)
- [Cardiac mechanics modeling](http://simvascular.github.io/docssvFSI-Structure.html)
- [Prescribed-motion LV modeling](http://simvascular.github.io/docssvFSI-Fluid.html#pres)

## Tutorial
- SimVascular group uploads hands-on tutorials to our [YouTube](https://www.youtube.com/channel/UCT61XgTRqpfb39Hyio9IqGQ) channel periodically. Here are some for `svFSI`:
  - Fluid-Structure Interaction (FSI): https://www.youtube.com/watch?v=QIpyThIAD7k
  - Cardiac electrophysiology modeling: https://www.youtube.com/watch?v=TCK3SmGwBa8
  - Cardiac mechanics modeling: https://www.youtube.com/watch?v=Jm3VSi6Aci8
- We also maintain a large collection of examples that showcase different functionalities of `svFSI`. You can find them here: https://github.com/SimVascular/svFSI-Tests. Each case includes a README file that explains the problem in hand and some key aspects of the software configuration.

## Pre/Post Processing Tool
We are also maintaining a collection of useful pre and post processing tools that are compatible with `svFSI`:
https://github.com/SimVascular/svFSI-Tools

## Contribute to `svFSI`
We welcome and appreciate all types of contributions to `svFSI`.
- Seek support, suggest new features or report bugs, please contact us through [GitHub Issues](https://github.com/SimVascular/svFSI/issues) or [SimTK forum](https://simtk.org/plugins/phpBB/indexPhpbb.php?f=188).
- Contribute your code to `svFSI`, please submit a pull request through GitHub.
- Share your novel applications of `svFSI` with the community, please consider contribute your case to [svFSI-Tests](https://github.com/SimVascular/svFSI-Tests).

## Citation
In preparation.
