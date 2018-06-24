#
# external location for installed vtk libs
#
# * can be built using scripts in "../Externals"
# * pre-built binaries for ubuntu 14.04 and gnu 4.8 compilers
#   can be downloaded using the "get-vtk-binaries.sh" script:
#    i.e. % ./get-vtk-binaries.sh ubuntu_14

OPEN_SOFTWARE_BINARIES_TOPLEVEL=${HOME}/svsolver/externals/bin/gnu-4.8/x64

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
