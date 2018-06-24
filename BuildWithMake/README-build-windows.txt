------------------------------------------------------------------------
            Compiling Instructions for svSolver on Windows
                       Revised 2016-09-18
------------------------------------------------------------------------

--------
Overview
--------

By default, svSolver is configured to build on linux using
makefiles.  You must override this to build on Windows.

Our base test configuration for linux is:

Windows 10
Intel 7 processor

Microsoft Visual Studio 2013 C++ compiler

and/or

ifort/icpc/icc intel 14 compilers

NOTE: The CMake system is currently broken for svSolver and
      is under revision.

-----------
Major Steps
-----------

1. Windows Prerequisities
-------------------------
The following packages are required to build svsolver

Cygwin 64-bit version (install in C:/cygwin64)
with required build tools (e.g. make, tclsh, zip)

Microsoft Visual Studio 2013 compiler

Intel Fortran Compiler

3. Checkout svSolver source code
--------------------------------
% git clone https://github.com/SimVascular/svSolver.git svsolver

2. vtk libraries
----------------
svSolver requires vtk libraries.  They can be build using the
the scripts in "../Externals", or pre-built binaries can be
downloaded using the script

% cd svsolver/BuildWithMake
% ./get-vtk-binaries.sh msvc-12.5

4. Override options
-------------------
Override defaults with:

  * cluster_overrides.mk
  * global_overrides.mk
  * site_overrides.mk
  * pkg_overrides.mk

See include.mk for all options.  Sample override files
can be found in:

SampleOverrides

to use one of these files, copy into local BuildWithMake
directory and modify as needed, e.g.:

cygwin64% cd svsolver/BuildWithMake
cygwin64% cp SampleOverrides/windows/global_overrides.mk .

6. Build
--------
cygwin64% cd svsolver/BuildWithMake
cygwin64% source CygwinHelpers/msvc_2013_x64
cygwin64% source CygwinHelpers/intel_fortan_14_x64
cygwin64% make

7. Running developer version
----------------------------
Binaries are in "BuildWithMake/Bin" directory.

8. Build release (NOTE: out-of-date!)
-----------------
cygwin64% cd svsolver/BuildWithMake/Release
cygwin64% make

9. Installing a distribution (NOTE: out-of-date!)
----------------------------
To be updated.
