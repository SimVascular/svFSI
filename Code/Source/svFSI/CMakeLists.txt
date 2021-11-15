# Copyright (c) Stanford University, The Regents of the University of
#               California, and others.
#
# All Rights Reserved.
#
# See Copyright-SimVascular.txt for additional details.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# add include path for headers on third party libraries
include_directories(${SV_SOURCE_DIR}/ThirdParty/parmetis_svfsi/simvascular_parmetis_svfsi/ParMETISLib)
include_directories(${SV_SOURCE_DIR}/ThirdParty/tetgen/simvascular_tetgen)
include_directories(${MPI_C_INCLUDE_PATH})

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-O3 -DNDEBUG -march=native")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -pthread -std=legacy")

  # Use below flags for debugging
  # set(CMAKE_Fortran_FLAGS "-O0")
  # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -pthread -std=legacy")
  # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -Wall -Wconversion -Wline-truncation -pedantic -fimplicit-none -fbacktrace -fbounds-check -p -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-O3 -march=native")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -pthread")
else()
  # nothing for now
  # may need to set for other compilers
endif()

# Find Trilinos package if requested
if(SV_USE_TRILINOS)

  # source: https://github.com/trilinos/Trilinos_tutorial/wiki/CMakeFindPackageTrilinosExample

  # It may be required to add a path for Trilinos
  # This can be accomplished by adding the Trilinos
  # -DCMAKE_PREFIX_PATH=~/sfw/Trilinos_build/lib/cmake/Trilinos
  FIND_PACKAGE(Trilinos)

  # If FIND_PACKAGE successfully found your Trilinos install, it will
  # set the Boolean flag Trilinos_FOUND.  The following IF statement
  # fails with a FATAL_ERROR if Trilinos was not found.  If it _was_
  # found, it prints out the values of some Trilinos configuration
  # details.  You may find them useful for building your application
  # that uses Trilinos.
  IF(Trilinos_FOUND)
    MESSAGE("\nFound Trilinos!  Here are the details: ")
    MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
    MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
    MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
    MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
    MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
    MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
    MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
    MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
    MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
    MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
    MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
    MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
    MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
    MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
    MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
    MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
    MESSAGE("End of Trilinos details\n")

    if(NOT Trilinos_LIBRARIES)
      MESSAGE(WARNING "Trilinos_FOUND is true but Trilinos_LIBRARIES is empty. This suggests that FIND_PACKAGE(Trilinos) failed quietly. Proceed with caution.")
    endif()

    MESSAGE("Setting WITH_TRILINOS to true\n")
    set(WITH_TRILINOS 1)
  ELSE()
    MESSAGE(WARNING "Could not find Trilinos. Compiling svFSI without Trilinos.")
  ENDIF()

endif()

# add trilinos flags and defines
if(WITH_TRILINOS)
  ADD_DEFINITIONS(-DWITH_TRILINOS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()


set(FSRCS TYPEMOD.f
         UTIL.f
         COMU.f
         CHNL.f
         LIST.f
         MATFUN.f
         CEPMOD_AP.f
         CEPMOD_FN.f
         CEPMOD_TTP.f
         CEPMOD_BO.f
         CEPMOD.f
         MOD.f
         ALLFUN.f
         BAFINI.f
         BF.f
         CEP.f
         CEPION.f
         CMM.f
         CONTACT.f
         DEBUG.f
         DISTRIBUTE.f
         EQASSEM.f
         FFT.f
         FLUID.f
         FS.f
         FSI.f
         HEATF.f
         HEATS.f
         IB.f
         INITIALIZE.f
         LELAS.f
         LHSA.f
         LOADMSH.f
         LOADNRB.f
         LS.f
         MAIN.f
         MATMODELS.f
         MESH.f
         NN.f
         NURBS.f
         OUTPUT.f
         PIC.f
         POST.f
         READFILES.f
         READMSH.f
         REMESH.f
         SETBC.f
         SHELLS.f
         STOKES.f
         STRUCT.f
         TXT.f
         USTRUCT.f
         VTKXML.f)

set(F90SRCS vtkComMod.f90
            vtkTypeParams.f90
            vtkXMLParser.f90
            vtkLegacyParser.f90)

set(CSRCS SPLIT.c
          vtkZpipe.c)

set(CXXSRCS remeshTet.cpp)


if(WITH_TRILINOS)
  set(CXXSRCS ${CXXSRCS} trilinos_linear_solver.cpp)

  # trilinos directories and libraries
  include_directories(${Trilinos_TPL_INCLUDE_DIRS})
  include_directories(${Trilinos_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS})
  link_directories(${Trilinos_TPL_LIBRARY_DIRS})
endif()


simvascular_add_executable(${SV_SVFSI_EXE}
  SRCS ${CSRCS}
  SRCS ${CXXSRCS}
  SRCS ${FSRCS}
  SRCS ${F90SRCS}
  DEV_SCRIPT_NAME "mysvfsi"
  INSTALL_SCRIPT_NAME "svfsi"
  INSTALL_COMP CoreExecutables
  INSTALL_DESTINATION ${SV_INSTALL_RUNTIME_DIR})


target_link_libraries(${SV_SVFSI_EXE}
  ${GLOBAL_LIBRARIES}
  ${INTELRUNTIME_LIBRARIES}
  ${ZLIB_LIBRARY}
  ${METIS_SVFSI_LIBRARY_NAME}
  ${PARMETIS_SVFSI_LIBRARY_NAME}
  ${TETGEN_LIBRARY_NAME}
  ${SV_MPI_Fortran_LIBRARIES}
  ${SV_LIB_SVFSILS_NAME}${SV_MPI_NAME_EXT}
  )

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # svFSI requires LAPACK
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  target_link_libraries(${SV_SVFSI_EXE} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_target_properties(svFSI PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(${SV_SVFSI_EXE} -mkl)
endif()  

# extra MPI libraries only if there are not set to NOT_FOUND or other null
if(SV_MPI_EXTRA_LIBRARY)
  target_link_libraries(${SV_SVFSI_EXE} ${SV_MPI_EXTRA_LIBRARY})
endif()

if(WITH_TRILINOS)
  target_link_libraries(${SV_SVFSI_EXE} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif()
