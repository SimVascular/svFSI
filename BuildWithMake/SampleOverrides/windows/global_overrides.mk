#
#  Override default cluster to build on Windows
#  with msvc 2013 compiler
# 

CLUSTER = x64_cygwin
CXX_COMPILER_VERSION = msvc-12.5
FORTRAN_COMPILER_VERSION = ifort

#
# external location for installed vtk libs
#
# * can be built using scripts in "../Externals"
# * pre-built binaries for Windows 10, vs2013 compilers
#   can be downloaded using the "get-vtk-binaries.sh" script:
#    i.e. % ./get-vtk-binaries.sh msvc-12.5

SVEXTERN_COMPILER_VERSION = $(CXX_COMPILER_VERSION)
OPEN_SOFTWARE_BINARIES_TOPLEVEL =C:/cygwin64/${HOME}/svsolver/externals/bin/$(SVEXTERN_COMPILER_VERSION)/x64
LICENSED_SOFTWARE_TOPLEVEL      =

#
# licensed is only needed if using leslib
#

LICENSED_SOFTWARE_TOPLEVEL= 

#
# Notes on MPI:
# * default is to use msmpi
# * can only build one at a time

SV_USE_DUMMY_MPI=0
