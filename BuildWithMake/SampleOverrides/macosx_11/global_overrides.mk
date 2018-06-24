#
#  Override default cluster to build on Mac OSX
#  with clang and intel fortran compiler
# 

CLUSTER = x64_macosx
CXX_COMPILER_VERSION = clang
FORTRAN_COMPILER_VERSION = ifort

#
# external location for installed vtk libs
#
# * can be built using scripts in "../Externals"
# * pre-built binaries for macosx clang 7.3 compilers
#   can be downloaded using the "get-vtk-binaries.sh" script:
#    i.e. % ./get-vtk-binaries.sh macosx_11

SVEXTERN_COMPILER_VERSION = clang-7.3
OPEN_SOFTWARE_BINARIES_TOPLEVEL=${HOME}/svsolver/externals/bin/$(SVEXTERN_COMPILER_VERSION)/x64

#
# licensed is only needed if using leslib
#

LICENSED_SOFTWARE_TOPLEVEL= 

#
# Notes on MPI:
# * default is to use dummy mpi
# * can only build one at a time
#

SV_USE_DUMMY_MPI=1
SV_USE_OPENMPI=0
SV_USE_MPICH=0
