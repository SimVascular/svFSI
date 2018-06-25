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
