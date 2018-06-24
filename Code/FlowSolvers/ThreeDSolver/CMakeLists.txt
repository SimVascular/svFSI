##=========================================================================
 #
 # Copyright (c) 2014-2015 The Regents of the University of California.
 # All Rights Reserved.
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
 #
 #=========================================================================#

set(INDENTS "${INDENTS}    ")

set(FLOWSOLVER_SUBDIRS ${FLOWSOLVER_SUBDIRS} svLS)

#-----------------------------------------------------------------------------
# Post Solver
#-----------------------------------------------------------------------------
if(SV_USE_SVPOST)
  set(FLOWSOLVER_SUBDIRS ${FLOWSOLVER_SUBDIRS} svPost)
  set(SV_USE_SOLVERIO ON)
endif()

#-----------------------------------------------------------------------------
# SolverIO
#-----------------------------------------------------------------------------
if(SV_USE_SOLVERIO)
  set(FLOWSOLVER_SUBDIRS ${FLOWSOLVER_SUBDIRS} SolverIO)
endif()

#-----------------------------------------------------------------------------
# Presolver
#-----------------------------------------------------------------------------
if(SV_USE_SVPRE)
  set(FLOWSOLVER_SUBDIRS ${FLOWSOLVER_SUBDIRS} svPre)
  set(SV_USE_SOLVERIO ON)
endif()

#-----------------------------------------------------------------------------
# Dummy MPI
#-----------------------------------------------------------------------------
if(SV_USE_DUMMY_MPI)
  set(FLOWSOLVER_SUBDIRS ${FLOWSOLVER_SUBDIRS} dummyMPI)
endif()

#-----------------------------------------------------------------------------
# svSolver
#-----------------------------------------------------------------------------
set(FLOWSOLVER_SUBDIRS ${FLOWSOLVER_SUBDIRS} svSolver)

#-----------------------------------------------------------------------------
# svFSI
#-----------------------------------------------------------------------------
if(NOT SV_USE_DUMMY_MPI)
  set(FLOWSOLVER_SUBDIRS ${FLOWSOLVER_SUBDIRS} svFSILS)
  set(FLOWSOLVER_SUBDIRS ${FLOWSOLVER_SUBDIRS} svFSI)
else()
  message(WARNING "svFSI requires MPI (not dummyMPI) and SV_USE_DUMMY_MPI \
is set to true. Not compiling svFSI or svFSILS target.")
endif()

#-----------------------------------------------------------------------------
# Optional header path
#-----------------------------------------------------------------------------
set(FLOWSOLVER_OPTIONS_HEADER_PATH ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "" FORCE)
dev_message("Flowsolver Header Path:  ${FLOWSOLVER_OPTIONS_HEADER_PATH}" )

#-----------------------------------------------------------------------------
# Solver feature support
#-----------------------------------------------------------------------------
include_directories("${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "${INDENTS}ThreeDSolver Options")
if((SV_USE_CORONARY EQUAL 0) AND (SV_USE_CLOSEDLOOP EQUAL 0) AND (VER_VARWALL EQUAL 0) AND (VER_USE_VTK EQUAL 0))
  message(STATUS "${INDENTS}    Compiling baseline version")
endif()
if(SV_USE_CORONARY)
  message(STATUS "${INDENTS}    Activating coronary bc feature")
endif()
if(SV_USE_CLOSEDLOOP)
  message(STATUS "${INDENTS}    Activating closed loop bc feature")
endif()
if(SV_USE_VARWALL)
  message(STATUS "${INDENTS}    Activating variable wall feature")
endif()
if(SV_USE_VTK)
  message(STATUS "${INDENTS}    Activating using VTK feature")
endif()

#-----------------------------------------------------------------------------
# Include and add subdirs
#-----------------------------------------------------------------------------
foreach(dir ${FLOWSOLVER_SUBDIRS})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${dir})
endforeach()

message(STATUS "${INDENTS}ThreeDSolver Components")
foreach(dir ${FLOWSOLVER_SUBDIRS})
  message(STATUS "${INDENTS}    ${dir}")
  add_subdirectory(${dir})
endforeach()

#-----------------------------------------------------------------------------
# Special option to also build a no mpi version!!
if(NOT SV_USE_DUMMY_MPI)
  if(SV_BUILD_ADDITIONAL_NOMPI_VERSION)
    set(SV_USE_DUMMY_MPI ON)
    set(SV_MPI_NAME_EXT "-nompi")
    include_directories(BEFORE dummyMPI)
    add_subdirectory(dummyMPI)
    add_subdirectory(svLS svLS${SV_MPI_NAME_EXT})
    add_subdirectory(svSolver svSolver${SV_MPI_NAME_EXT})
    set(SV_USE_DUMMY_MPI OFF)
  endif()
endif()
#-----------------------------------------------------------------------------


