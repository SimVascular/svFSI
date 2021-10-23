
This file describes a general procedure to install `svFSI` dependencies. While most dependencies provided below are built from source, users are encouraged to use package managers such as `apt` for Ubuntu and `brew` for Mac OSX when available. Users may modify the parameters according to the local module and compiler environment. Please note that the below procedures are primarily intended for moderate to advanced Linux users only. Others are encouraged to use precompiled binaries at [SimTK](https://simtk.org/frs/index.php?group_id=188) or the `Quick Build` procedure outlined in the README.md file.

<hr style="border:2px solid gray"> </hr>

## Intel oneAPI Toolkits

Intel has release its next-gen software development suite, oneAPI, free of charge. This suite provides comprehensive compiler support for Linux, Windows and limited support for macOS. To build `svFSI` with Intel oneAPI compilers, the following two toolkits should be installed:
[Intel oneAPI Base Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html#gs.cld9rl) and [Intel oneAPI HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit.html#gs.cldcfd). Detailed installation instructions are provided on the download page of each toolkit.

Intel oneAPI Base Toolkit provides development tools for general computing and consists of over a dozen components. Users only need to install selected few in order to build `svFSI`. We recommend the following:
- Intel oneAPI Collective Communications Library
- Intel oneAPI DPC++/C++ Compiler
- Intel oneAPI DPC++ Library
- Intel oneAPI Math Kernel Library
- Intel oneAPI Threading Building Blocks
- Intel Distribution for GDB
- Intel DPC++ Compatibility Tool
- Intel VTune Profiler

Intel oneAPI HPC Toolkit complements the Base Toolkit and provides necessary components for high-performance computing, such as Fortran compiler and MPI library. We recommend install all of the components. Please note that there is no Intel MPI library for macOS as of version 2021.3.0.

<hr style="border:2px solid gray"> </hr>

## Curses CMake

Before installing cmake, install libncurses5-dev from apt as,

```
sudo apt-get install libncurses5-dev
```

to get ccmake installed

<hr style="border:2px solid gray"> </hr>

## CMAKE 3.20.5

Downloaded cmake version 3.20.5 from,

https://github.com/Kitware/CMake/releases/tag/v3.20.5

Build using the following commands:

```bash
sudo mkdir -p /opt/cmake/3.20.5
tar -zxvf cmake-3.20.5.tar.gz
cd cmake-3.20.5
mkdir build
cd build/
../bootstrap  --prefix=/opt/cmake/3.20.5
make
sudo make install
```

<hr style="border:2px solid gray"> </hr>

## MPICH

Downloaded mpich version 3.4.2 from,

https://www.mpich.org/downloads/

Build using the following commands:

```bash
sudo mkdir /opt/mpich/3.4.2
tar -zxvf mpich-3.4.2.tar.gz
mkdir build
cd build/
../mpich-3.4.2/configure  --prefix=/opt/mpich/3.4.2 --enable-fast=O3,ndebug --without-timing  --without-mpit-pvars  --with-device=ch3:nemesis
make
sudo make install
```

<hr style="border:2px solid gray"> </hr>

## OpenMPI

Downloaded openmpi version 4.1.1 from,

https://www.open-mpi.org/software/ompi/v4.1/

Build using the following commands:

```bash
sudo mkdir -p /opt/openmpi/4.1.1
tar -zxvf openmpi-4.1.1.tar.gz
mkdir build
cd build/
../openmpi-4.1.1/configure --prefix=/opt/openmpi/4.1.1 --enable-static
make
sudo make install
```

<hr style="border:2px solid gray"> </hr>

## LAPACK

Downloaded lapack version 3.9.1 from,

http://www.netlib.org/lapack/#_lapack_version_3_9_1

### Build shared libraries using the following commands:

```bash
sudo mkdir -p /opt/netlib/lapack/3.9.1/shared
tar -zxvf lapack-3.9.1.tar.gz
mkdir build-shared
cd build-shared/
ccmake -DCMAKE_PREFIX_PATH:PATH=/opt/netlib/lapack/3.9.1/shared \
   -DCMAKE_BUILD_TYPE:STRING="Release" \
   -DCMAKE_C_FLAGS:STRING="-O3 -DNDEBUG -march=native" \
   -DCMAKE_Fortran_FLAGS:STRING="-O3 -DNDEBUG -march=native -frecursive" \
   -DBUILD_SHARED_LIBS:BOOL=ON \
   -DBUILD_DEPRECATED:BOOL=ON \
   ../lapack-3.9.1
make
sudo make install
```

Create a link to the library at a standard location such as

```bash
sudo ln -sf /opt/netlib/lapack/3.9.1/static/liblapack.so  /usr/local/lib
sudo ln -sf /opt/netlib/lapack/3.9.1/static/libblas.so  /usr/local/lib
```

### Build static libraries using the following commands:

```bash
sudo mkdir -p /opt/netlib/lapack/3.9.1/static
tar -zxvf lapack-3.9.1.tar.gz
mkdir build-static
cd build-static/
ccmake -DCMAKE_PREFIX_PATH:PATH=/opt/netlib/lapack/3.9.1/static \
   -DCMAKE_BUILD_TYPE:STRING="Release" \
   -DCMAKE_C_FLAGS:STRING="-O3 -DNDEBUG -march=native" \
   -DCMAKE_Fortran_FLAGS:STRING="-O3 -DNDEBUG -march=native -frecursive" \
   -DBUILD_DEPRECATED:BOOL=ON \
   ../lapack-3.9.1
make
sudo make install
```

Create a link to the library at standard location

```bash
sudo ln -sf /opt/netlib/lapack/3.9.1/static/liblapack.a  /usr/local/lib
sudo ln -sf /opt/netlib/lapack/3.9.1/static/libblas.a  /usr/local/lib
```

<hr style="border:2px solid gray"> </hr>

## Boost

Boost installed version: boost_1_66_0.tar.gz

Build using the following commands:

```bash
sudo mkdir /opt/boost/1.66.0
tar -zxvf boost_1_66_0.tar.gz
cd boost_1_66_0
./bootstrap.sh --prefix=/opt/boost/1.66.0 --with-libraries=all
sudo ./b2 install
```

<hr style="border:2px solid gray"> </hr>

## HDF5

HDF5 installed version: hdf5-1.10.4.tar.gz

Build using the following commands:

```bash
sudo mkdir -p /opt/hdf5/1.10.4
tar -zxvf hdf5-1.10.4.tar.gz
cd hdf5-1.10.4
mkdir build
cd build

make
sudo make install
```

<hr style="border:2px solid gray"> </hr>

## HYPRE

HYPRE installed version: hypre-2.22.0.tar.gz
Downloaded from https://github.com/hypre-space/hypre/tags

Built shared libraries using the following commands:

```bash
sudo mkdir -p /opt/hypre/2.22.0/gnu-mpich
tar -zxvf hypre-2.22.0.tar.gz
cd hypre-2.22.0/src
./configure --prefix=/opt/hypre/2.22.0/gnu-mpich/shared --enable-shared  --with-extra-CFLAGS=" -march=native" --with-extra-CXXFLAGS=" -march=native" --with-MPI --with-lapack --with-blas
sudo make install
```

<hr style="border:2px solid gray"> </hr>

## MUMPS

MUMPS installed version: MUMPS_5.4.0.tar.gz
Downloaded via http://mumps.enseeiht.fr/index.php?page=dwnld

### Unpack

Unpack the tar ball and set folder structure:

```bash
sudo mkdir -p /opt/mumps/5.4.0
tar -zxvf MUMPS_5.4.0.tar.gz
cd MUMPS_5.4.0
```

### Prepare Makefile.inc

Copy a template Makefile.inc to the src/ directory as,

```bash
cp Make.inc/Makefile.debian.PAR ./Makefile.inc
```

Modify the variables in Makefile.inc according to the local environment, libraries, etc. For e.g., we use the following settings for parallel a build of MUMPS with metis, scotch, lapack, blas, scalapack, and mpich:

```bash
#
#  This file is part of MUMPS 5.4.0, released
#  on Tue Apr 13 15:26:30 UTC 2021
#
# These settings for an Ubuntu PC with custom-built packages :
# metis (parmetis), scotch (ptscotch), lapack, blas, scalapack, mpich, gfortran

# Begin orderings
SCOTCHDIR  = /opt/scotch/6.1.0
ISCOTCH    = -I/opt/scotch/6.1.0/include
LSCOTCH    = -L$(SCOTCHDIR)/lib -lptesmumps -lptscotch -lptscotcherr -lscotch

PORDDIR    = $(topdir)/PORD
IPORD      = -I$(PORDDIR)/include
LPORDDIR   = $(PORDDIR)/lib
LPORD      = -L$(LPORDDIR) -lpord

METISDIR   = /opt/metis/5.1.0
PMETISDIR  = /opt/parmetis/4.0.3
IMETIS     = -I$(PMETISDIR)/include -I$(METISDIR)/include
LMETIS     = -L$(PMETISDIR)/lib -lparmetis -L$(METISDIR)/lib -lmetis

# Corresponding variables reused later
ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis
ORDERINGSC = $(ORDERINGSF)

LORDERINGS  = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)
# End orderings
################################################################################

PLAT    =
LIBEXT  = .a
OUTC    = -o
OUTF    = -o
RM      = /bin/rm -f
CC      = mpicc
FC      = mpif90
FL      = mpif90
AR      = ar vr
RANLIB  = ranlib
LAPACK  = /opt/netlib/lapack/3.9.1/static/lib/liblapack.a
SCALAP  = /opt/netlib/scalapack/2.1.0/libscalapack.a

INCPAR = -I/opt/mpich/3.4.2/include

LIBPAR = $(SCALAP) $(LAPACK) -L/opt/mpich/3.4.2/lib -lfmpich -lmpich -lmpi

INCSEQ = -I$(topdir)/libseq
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq -lmpiseq

LIBBLAS = /opt/netlib/lapack/3.9.1/static/lib/libblas.a
LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
OPTF    = -O3 -DNDEBUG -march=native -fopenmp
OPTL    = -O3 -DNDEBUG -march=native -fopenmp
OPTC    = -O3 -DNDEBUG -march=native -fopenmp
#End Optimized options

INCS = $(INCPAR)
LIBS = $(LIBPAR)
LIBSEQNEEDED =
```

### Compile

```bash
make all
```

### Build shared libraries to link with Trilinos

To build shared libraries of MUMPS, the following steps may help:

- Add the "-fPIC" option to your compiler options in Makefile.inc
- Run `make clean` and recompile all MUMPS source files, to build the MUMPS ".a" libraries
- Run something like `ld -shared -o libdmumps.so libdmumps.a` to create a shared library "libdmumps.so".
- Repeat the last step for all the static libraries in the lib/ folder to create the corresponding shared libraries.


### Copy to the destination folder

Copy `include`, `lib`, and `examples` to the destination folder.

```bash
sudo cp -r include/ lib/ examples/ /opt/mumps/5.4.0
```

<hr style="border:2px solid gray"> </hr>

## Trilinos

Trilinos installed version: Trilinos-trilinos-release-13-0-1.tar.gz

Downloaded from
https://github.com/trilinos/Trilinos/releases/tag/trilinos-release-13-0-1

Trilinos preconfiguration steps include:

```bash
sudo mkdir -p /opt/trilinos/13.0.1/gnu-mpich
tar -zxvf Trilinos-trilinos-release-13-0-1.tar.gz
mkdir build_mpich
cd build_mpich
```

Compiler flags for performance

```bash
CXXFLAGS = -O3 -DNDEBUG -march=native
CFLAGS = -O3 -DNDEBUG -march=native
FCFLAGS = -O3 -march=native
```

CMake configuration for GNU compilers:

```bash
    cmake                                            \
    -DTrilinos_ENABLE_Amesos=ON                      \
    -DTrilinos_ENABLE_AztecOO=ON                     \
    -DTrilinos_ENABLE_Epetra=ON                      \
    -DTrilinos_ENABLE_EpetraExt=ON                   \
    -DTrilinos_ENABLE_Ifpack=ON                      \
    -DTrilinos_ENABLE_ML=ON                          \
    -DTrilinos_ENABLE_MueLu=ON                       \
    -DTrilinos_ENABLE_ROL=ON                         \
    -DTrilinos_ENABLE_Sacado=ON                      \
    -DTrilinos_ENABLE_Teuchos=ON                     \
    -DTrilinos_ENABLE_Zoltan=ON                      \
    -DTrilinos_VERBOSE_CONFIGURE=OFF                 \
    -DTPL_ENABLE_Boost:BOOL=ON                       \
    -DTPL_ENABLE_BLAS:BOOL=ON                        \
    -DTPL_ENABLE_HDF5:BOOL=ON                        \
        -DHDF5_INCLUDE_DIRS:PATH=<PATH TO HDF5 INSTALL INCLUDE DIR> \
        -DHDF5_LIBRARY_DIRS:PATH=<PATH TO HDF5 INSTALL LIB DIR>     \
    -DTPL_ENABLE_HYPRE:BOOL=ON                       \
    -DTPL_ENABLE_LAPACK:BOOL=ON                      \
    -DTPL_ENABLE_MPI=ON                              \
        -DMPI_USE_COMPILER_WRAPPERS=ON                        \
        -DMPI_C_COMPILER:PATH="${MPI_BIN}/mpicc"              \
        -DMPI_CXX_COMPILER:PATH="${MPI_BIN}/mpic++"           \
        -DMPI_Fortran_COMPILER:PATH="${MPI_BIN}/mpif77"       \
        -DCMAKE_C_COMPILER:PATH="${MPI_BIN}/mpicc"            \
        -DCMAKE_CXX_COMPILER:PATH="${MPI_BIN}/mpic++"         \
        -DCMAKE_Fortran_COMPILER:PATH="${MPI_BIN}/mpif77"     \
        -DMPI_BASE_DIR="${MPI_DIR}"                  \
    -DTPL_ENABLE_MUMPS:BOOL=ON                       \
        -DMUMPS_INCLUDE_DIRS:PATH=<PATH TO MUMPS INSTALL INCLUDE DIR> \
        -DMUMPS_LIBRARY_DIRS:PATH=<PATH TO MUMPS INSTALL LIB DIR>    \
    -DBUILD_SHARED_LIBS=ON                           \
    -DCMAKE_VERBOSE_MAKEFILE=OFF                     \
    -DCMAKE_BUILD_TYPE=RELEASE                       \
    -DCMAKE_INSTALL_PREFIX:PATH=<PATH TO TRILINOS INSTALL DIR> \
    ../
```

Minimalist CMake configuration for Intel oneAPI Toolkits
```bash
#!/bin/bash

INTEL_MPI_DIR=/opt/pkg/intel/oneapi/mpi/2021.3.0
INTEL_MKL_DIR=/opt/pkg/intel/oneapi/mkl/2021.3.0

cmake \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/pkg/trilinos/13.0.1_intel \
    -D TPL_ENABLE_MPI=ON \
    -D MPI_BASE_DIR="${INTEL_MPI_DIR}" \
    -D TPL_ENABLE_MKL:BOOL=ON \
    -D MKL_LIBRARY_DIRS:FILEPATH="${INTEL_MKL_DIR}/lib/intel64" \
    -D MKL_INCLUDE_DIRS:FILEPATH="${INTEL_MKL_DIR}/include" \
    -D Trilinos_ENABLE_OpenMP=OFF \
    -DTPL_ENABLE_BLAS:BOOL=ON \
    -D BLAS_LIBRARY_DIRS:FILEPATH="${INTEL_MKL_DIR}/lib/intel64" \
    -D BLAS_LIBRARY_NAMES:STRING="mkl_intel_lp64;mkl_sequential;mkl_core" \
    -DTPL_ENABLE_LAPACK:BOOL=ON \
    -D LAPACK_LIBRARY_DIRS:FILEPATH="${INTEL_MKL_DIR}/lib/intel64" \
    -D LAPACK_LIBRARY_NAMES:STRING="mkl_intel_lp64;mkl_sequential;mkl_core" \
    -DTrilinos_ENABLE_FLOAT=ON \
    -DTrilinos_ENABLE_TESTS:BOOL=ON \
    -DTrilinos_ENABLE_Epetra:BOOL=ON \
    -DTrilinos_ENABLE_AztecOO:BOOL=ON \
    -DTrilinos_ENABLE_Amesos:BOOL= ON \
    -DTrilinos_ENABLE_ML:Bool=ON \
    -DTrilinos_ENABLE_MueLu:Bool=ON \
    -DTrilinos_ENABLE_Ifpack:BOOL=ON \
    -DCMAKE_BUILD_TYPE=RELEASE \
    -D BUILD_SHARED_LIBS=ON \
    ../
```
Here, `BLAS_LIBRARY_NAMES` and `LAPACK_LIBRARY_NAMES` are platform/architecture specific. Please consult [Intel Link Line Advisor](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html).

Finally, compile and install

```bash
make -j <num_processes>
sudo make install
```

<hr style="border:2px solid gray"> </hr>


