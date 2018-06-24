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

# Set c++ sources
set(CXXSRCS cmd.cxx displacements.cxx helpers.cxx supre-cmds.cxx supre.cxx
  presolver-vtk-cmds.cxx)

# Add c++ sources for variable wall module
if(SV_USE_VARWALL)
  set(CXXSRCS ${CXXSRCS} Laplace_Solve.cxx Transient_Laplace_Solve.cxx displacements_var_prop.cxx)
endif()

# Add c sources
set(HDRS cmd.h)
set(CSRCS directsolve.c)

# Add fortran sources
set(FSRCS stannspcg.f)

# add fortran sources for variable wall module
if(SV_USE_VARWALL)
  set(FSRCS ${FSRCS} gmresfortran.f90)
endif()

add_library(${SV_LIB_PRESOLVER_NAME} ${SV_LIBRARY_TYPE} ${FSRCS})

simvascular_add_executable(${SV_PRESOLVER_EXE}
  SRCS ${CXXSRCS} ${CSRCS}
  DEV_SCRIPT_NAME "mypre"
  INSTALL_SCRIPT_NAME "svpre"
  INSTALL_COMP CoreExecutables
  INSTALL_DESTINATION ${SV_INSTALL_RUNTIME_DIR})

target_link_libraries(${SV_PRESOLVER_EXE} ${SV_LIB_PRESOLVER_NAME}
  ${GLOBAL_LIBRARIES} ${INTELRUNTIME_LIBRARIES}
  ${METIS_LIBRARY} ${SPARSE_LIBRARY}
  ${NSPCG_LIBRARY} ${VTK_LIBRARIES}
  ${SV_LIB_SOLVERIO_NAME})

target_link_libraries(${SV_LIB_PRESOLVER_NAME} ${NSPCG_LIBRARY})

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set_target_properties(${SV_PRESOLVER_EXE} PROPERTIES LINK_FLAGS -nofor-main)
endif()

if(SV_INSTALL_LIBS)
  install(TARGETS ${SV_LIB_PRESOLVER_NAME}
    RUNTIME DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT CoreExecutables
    LIBRARY DESTINATION ${SV_INSTALL_LIBRARY_DIR} COMPONENT CoreLibraries
    ARCHIVE DESTINATION ${SV_INSTALL_ARCHIVE_DIR} COMPONENT CoreLibraries
    )
endif()
if(SV_INSTALL_HEADERS)
  install(FILES ${HDRS}
    DESTINATION ${SV_INSTALL_INCLUDE_DIR}/core COMPONENT CoreHeaders
    )
endif()


