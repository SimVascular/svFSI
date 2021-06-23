
This file describes a generaly procedure to install `svFSI` dependencies. While most dependencies provided below are built from source, users are encouraged to use package managers such as `apt` for Ubuntu and `brew` for Mac OSX when available. Users may modify the parameters according to the local module and compiler environment. Please note that the below procedures are primarily intended for moderate to advanced linux users only. Others are encouraged to use precompiled binaries at [SimTK](https://simtk.org/frs/index.php?group_id=188) or the `Quick Build` procedure outlined in the README.md file.

## =================================================================
## Curses CMake

Before installing cmake, install libncurses5-dev from apt as,

```
sudo apt-get install libncurses5-dev
```

to get ccmake installed

## =================================================================
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

## =================================================================
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

## =================================================================
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

# =================================================================
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

Create a link to the library at standard location

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

## =================================================================
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

# =================================================================
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

## =================================================================
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

## =================================================================
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

### Compiler flags for performance

CXXFLAGS = -O3 -DNDEBUG -march=native
CFLAGS = -O3 -DNDEBUG -march=native
FCFLAGS = -O3 -march=native

### CMake Configuration:

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
    -DTPL_ENABLE_HDF5:BOOL=ON                       \
        -DHDF5_INCLUDE_DIRS:PATH=<PATH TO HDF5 INSTALL INCLUDE DIR> \
        -DHDF5_LIBRARY_DIRS:PATH=<PATH TO HDF5 INSTALL LIB DIR>    \
    -DTPL_ENABLE_HYPRE:BOOL=ON                       \
    -DTPL_ENABLE_LAPACK:BOOL=ON                      \
    -DTPL_ENABLE_MPI=ON                              \
        -DMPI_USE_COMPILER_WRAPPERS=ON                  \
        -DMPI_C_COMPILER:PATH="${MPI_BIN}/mpicc"          \
        -DMPI_CXX_COMPILER:PATH="${MPI_BIN}/mpic++"           \
        -DMPI_Fortran_COMPILER:PATH="${MPI_BIN}/mpif77"       \
        -DCMAKE_C_COMPILER:PATH="${MPI_BIN}/mpicc"            \
        -DCMAKE_CXX_COMPILER:PATH="${MPI_BIN}/mpic++"         \
        -DCMAKE_Fortran_COMPILER:PATH="${MPI_BIN}/mpif77"     \
        -DMPI_BASE_DIR="${MPI_DIR}"
    -DBUILD_SHARED_LIBS=ON                           \
    -DCMAKE_VERBOSE_MAKEFILE=OFF                     \
    -DCMAKE_BUILD_TYPE=RELEASE                       \
    -DCMAKE_INSTALL_PREFIX:PATH=<PATH TO TRILINOS INSTALL DIR> \
    ../
```

Finally, compile and install

```bash
make -j 4
sudo make install
```

## =================================================================

