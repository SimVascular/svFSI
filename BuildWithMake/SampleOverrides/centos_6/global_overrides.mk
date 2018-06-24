#
# external location for installed vtk libs
#
# * can be built using scripts in "../Externals"
# * pre-built binaries for centos 6.8, gnu 4.4.7 compilers
#   can be downloaded using the "get-vtk-binaries.sh" script:
#    i.e. % ./get-vtk-binaries.sh centos_6

OPEN_SOFTWARE_BINARIES_TOPLEVEL=${HOME}/svsolver/externals/bin/gnu-4.4/x64

#
# licensed is only needed if using leslib
#

LICENSED_SOFTWARE_TOPLEVEL= 

#
# Notes on MPI:
# * default is to use dummy mpi
# * can only build one at a time
# * use modules on centos: 
#    e.g.  % module add openmpi-x86_64
#    e.g.  % module add mpich-x864_64
#

SV_USE_DUMMY_MPI=1
SV_USE_OPENMPI=0
SV_USE_MPICH=0
