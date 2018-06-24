------------------------------------------------------------------------
            Compiling Instructions for svSolver on Mac OSX
                       Revised 2016-09-18
------------------------------------------------------------------------

--------
Overview
--------

By default, svSolver is configured to build on osx using
makefiles.  Our base test configuration for OSX is:

Apple OSX 10.11 64-bit
Intel 7 processor

clang/clang++ version 7.3
macports mpich/gfortran5

and/or

ifort/icpc/icc intel compilers

NOTE: The CMake system is currently broken for svSolver and
      is under revision.

-----------
Major Steps
-----------

1. System Prerequisities
------------------------
The following packages are required to build svsolver

XCode command line tools are required
% xcode-select --install
% sudo xcodebuild -license

MacPorts is required which can be downloaded at https://www.macports.org
The following packages are required to build simvascular

### compilers
% sudo port install gcc5
% sudo port install gfortran5
% sudo port install mpich-devel-clang

2. Checkout svSolver source code
--------------------------------
% git clone https://github.com/SimVascular/svSolver.git svsolver

3. vtk libraries
----------------
svSolver requires vtk libraries.  They can be build using the
the scripts in "../Externals", or pre-built binaries can be
downloaded using the script

% cd svsolver/BuildWithMake
% ./get-vtk-binaries.sh macosx_11

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

% cd svsolver/BuildWithMake
% cp SampleOverrides/macosx_11/global_overrides.mk

6. Build
--------
% cd svsolver/BuildWithMake
% make

7. Running developer version
----------------------------
Binaries are in "BuildWithMake/Bin" directory.

You may need to set the path for the intel shared libraries before running svsolver, e.g.:

% export DYLD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib/

8. Build release (NOTE: out-of-date!)
-----------------
To be updated.

9. Installing a distribution (NOTE: out-of-date!)
----------------------------
To be updated.
