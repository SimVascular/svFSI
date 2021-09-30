
# **Build svFSI from Source**

Below are the instructions to build svFSI on Ubuntu and MacOS.

<hr style="border:2px solid gray"> </hr>

## svFSI on Ubuntu 18.04

### Download svFSI

`svFSI` source can be downloaded from GitHub at:
https://github.com/SimVascular/svFSI

### Dependencies
Dependencies for `svFSI` include:

- C, Fortran compilers (gcc/gfortran, icc/ifort,..)
- MPI (mpich, openmpi, intel MPI,..)
- CMake
- LAPACK and BLAS
- Trilinos (optional)

Please refer to [`INSTALL-DEPS.md`](./INSTALL-DEPS.md) for more information on installing dependencies.

Recommended folder structure:

```bash
mkdir svFSI
cd svFSI
git clone https://github.com/SimVascular/svFSI.git
mv svFSI src
mkdir build
cd build
```

If you have `ssh` keys set up on GitHub, you may clone using the following command instead of the above `https` ones,

```bash
git clone git@github.com:SimVascular/svFSI.git
```

This structure creates a separate directory for `svFSI` git repository `src` and a separate `build` folder where `svFSI` is compiled.

### Build

Tested environment: Ubuntu 18.04,
gcc/g++/gfortran 7.5.0,
mpich 3.4.0,
cmake 3.20.5,
lapack/blas 3.9.1
   ```bash
   ccmake ..
   ```

Or

Tested environment: Ubuntu 18.04, Intel oneAPI Base Toolkit 2021.3.0, Intel oneAPI HPC Toolkit 2021.3.0
   ```bash
   ccmake \
​      -DCMAKE_C_COMPILER:PATH=${ONEAPI_ROOT}/compiler/2021.3.0/linux/bin/intel64/icc \
​      -DCMAKE_CXX_COMPILER:PATH=${ONEAPI_ROOT}/compiler/2021.3.0/linux/bin/intel64/icpc \
​      -DCMAKE_Fortran_COMPILER:PATH=${ONEAPI_ROOT}/compiler/2021.3.0/linux/bin/intel64/ifort ..
   ```
   Here `ONEAPI_ROOT` is an environment variable automatically set by Intel oneAPI.

### Build with Trilinos (optional)

Tested environment:
Trilinos & its dependencies(
    boost 1.66.0,
    hdf5 1.10.4,
    hypre 2.22.0,
    trilinos 13.01)

To compile `svFSI` with Trilinos, the following changes need to be made to the file, `Code/CMake/SimVascularOptions.cmake`

```bash
option(SV_USE_TRILINOS "Use Trilinos Library with svFSI" ON)
```

Further, a `CMAKE_PREFIX_PATH` should be provided as command line argument pointing to the Trilinos library as,

```bash
ccmake -DCMAKE_PREFIX_PATH:PATH=$(TRILINOS_DIR)/lib/cmake/Trilinos <path_to_svFSI_source>
```

where `TRILINOS_DIR` should be replaced by the path to the Trilinos installation directory such as `/opt/trilinos/13.01/gnu-mpich`

### Compiler flags for performance

To build `svFSI` for performance, the following CMake options may be used,

```bash
CMAKE_BUILD_TYPE:STRING = "RELEASE"
SV_BUILD_TYPE:STRING = "REEASE"
CMAKE_C_FLAGS:STRING = "-O3 -DNDEBUG -march=native"
CMAKE_CXX_FLAGS:STRING = "-O3 -DNDEBUG -march=native"
```

The CMake Fortran flags for `svFSI` are hard-coded. So it is the user's responsibility to use optimized flags depending on the compiler environment. These changes need to be made in the source code as,

open file `Code/Source/svFSI/CMakeLists.txt` and modify `CMAKE_Fortran_FLAGS` as,

```bash
  set(CMAKE_Fortran_FLAGS "-O3 -DNDEBUG -march=native")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -pthread -std=legacy")
```

and for the linear solver `svFSILS`, modify `CMAKE_Fortran_FLAGS` in the file `Code/Source/svFSILS/CMakeLists.txt` as,

```bash
  set(CMAKE_Fortran_FLAGS "-O3 -DNDEBUG -march=native")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -pthread -std=legacy")
```

Note that the option `std=legacy` is used for more recent versions of gcc (version >= 10).

### Additional CMake settings

Occasionally, C and CXX compilers may be wrongly identified by the CMake system used for `svFSI`. By default, the C compiler used is `/usr/bin/cc` and the CXX compiler is `/usr/bin/c++`. However, depending on the modules loaded on HPC or one's local compiler environment, these may not be pointing to the desired compilers such as those wrapped by mpicc or mpicxx.

In such situations, we recommend providing the correct C and CXX compilers using the CMake variables as,

```bash
ccmake -DCMAKE_C_COMPILER:PATH=/usr/bin/gcc -DCMAKE_CXX_COMPILER:PATH=/usr/bin/g++  <path_to_svFSI_source>
```

### CMake configuration:

An example CMake command for compiling `svFSI` with Trilinos is given below:

```bash
    cmake                                                           \
    -DCMAKE_C_COMPILER:PATH=/usr/bin/gcc                            \
    -DCMAKE_CXX_COMPILER:PATH=/usr/bin/g++                          \
    -DCMAKE_PREFIX_PATH:PATH=$(TRILINOS_DIR)/lib/cmake/Trilinos     \
    -DCMAKE_BUILD_TYPE:STRING="Release"                             \
    -DCMAKE_C_FLAGS:STRING="-O3 -DNDEBUG -march=native"             \
    -DCMAKE_CXX_FLAGS:STRING="-O3 -DNDEBUG -march=native"           \
    -DSV_BUILD_TYPE_DIR:STRING="Release"             \
    <../svFSI_src>
```

Note that the `-march` flag could be chosen depending on the processor. A variety of options are available for Intel and AMD processors at,

https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html

Finally, run `make` for compiling `svFSI`.

If the compilation proceeds successfully, a binary file is created in `svFSI-build/bin/svFSI`.

You may quickly run an `ldd` to check if all the required libraries are available such as,

```bash
ldd svFSI-build/bin/svFSI
```
The code will not run if any library is not found. A wrapper script is also available in `svFSI-build/mysvFSI` which sets environmental variables.

<hr style="border:2px solid gray"> </hr>

## `svFSI` on Mac OSX

Follow the steps below to install `svFSI` on Mac.

1. Install Xcode command-line tools using the command below:

```bash
xcode-select install
```

2. Check if the default terminal is csh or bash. You may open the Terminal/Preferences and change the default shell to bash (/bin/bash). This is recommended as we may add some environmental variables later.

3. Install [Homebrew](https://docs.brew.sh/Installation)

4. Check if `/usr/local/bin` appears before `/usr/bin` in the `PATH` environmental variable.

```bash
echo $PATH
```

If not, open `~/.bash_profile` and add the following line:

```bash
export PATH="/usr/local/bin:$PATH"
```

5. Install the following packages using Homebrew

```bash
brew install gcc
brew install mpich
brew install cmake
```

Verify if MPI compilers are pointing to the correct compilers as

```bash
which mpicc
which mpicxx
which mpif90
```

The CC and CXX compilers should point to the default `clang` compilers, while the FC/F90 compiler should point to `gfortran`. Set the environmental variable `LIBRARY_PATH` to the location of MPI_LIBRARY (usually `/usr/local/lib`) as,

```bash
export LIBRARY_PATH="/usr/local/lib:$LIBRARY_PATH"
```

in `~/.bash_profile`.

6. Mac OS usually comes with default `blas` and `lapack` packages. If you prefer to use `OpenBLAS`, you may install it using Homebrew as,

```bash
brew install openblas
brew install lapack
```

You should also set the below environmental variables in `~/.bash_profile` as,

```bash
BLASDIR="/usr/local/Cellar/openblas/0.3.15_1"
export LIBRARY_PATH="$BLASDIR/lib:$LIBRARY_PATH"
```

Please make sure to update the `BLASDIR` location depending on the version of OpenBLAS you install.

7. Download `svFSI` from GitHub:

```bash
mkdir svFSI
cd svFSI
git clone https://github.com/SimVascular/svFSI.git
mv svFSI/  src
```

If you have `ssh` keys set up on GitHub, you may clone using the following command instead of the above `https` ones,

```bash
git clone git@github.com:SimVascular/svFSI.git
```

8. Create a separate `build` folder and compile using the following commands:

```bash
mkdir build
cd build
ccmake ../src
make
```
<hr style="border:2px solid gray"> </hr>