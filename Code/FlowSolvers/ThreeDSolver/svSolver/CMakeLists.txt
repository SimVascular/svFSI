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

if(MSVC)
  string(REGEX REPLACE "STACK:10000000,10000000" "STACK:512000000,512000000" CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
endif()
# Include c++ files
set(CXXSRCS Cinput.cxx ValType.cxx input_fform.cxx partition.cxx main.cxx)

# Include c files
set(CSRCS getIntPnts.c ShapeFunctions_uniformP.c shptet.c symtri.c symtet.c) # for shape functions
set(CSRCS ${CSRCS} write_restart_files.c) # solver io files
set(CSRCS ${CSRCS} tmrc.c) # this is dummy timing function and should go away

if(NOT (LINUX AND IS64))
  set(CSRCS ${CSRCS} missing_win32.c)
endif()

set(FSRCS modules.f)

if(SV_USE_CORONARY)
  set(FSRCS ${FSRCS} modules-coronary.f)
endif()

if (SV_USE_CLOSEDLOOP)
  set(FSRCS ${FSRCS} modules-closedloop.f)
endif()

set(FSRCS ${FSRCS} solver_subroutines.f)

set(FSRCS ${FSRCS} addelmpvsqforsvls-closedloop-coronary.f
                   asbmfg-closedloop-varwall.f
                   asbnabi-closedloop.f
                   asbnasc-closedloop.f
                   asigmr-closedloop.f
                   elmgmr-closedloop-varwall.f
                   elmpvsq-closedloop-coronary.f
                   genbkb-varwall.f
                   genini-coronary.f
                   itrdrv-closedloop-varwall-coronary.f
                   proces-varwall.f
                   readnblk-varwall.f
                   updhistconv-coronary.f)

if(SV_USE_CORONARY)
  set(FSRCS ${FSRCS} coronary_subroutines.f)
endif()

# Include files for the CLOSEDLOOP VERSION
if(SV_USE_CLOSEDLOOP)
  set(FSRCS ${FSRCS} closedloop_subroutines.f)
endif()

# Include files for the VARWALL VERSION
if (SV_USE_VARWALL)
  set(FSRCS ${FSRCS} varwall_subroutines.f)
endif()

# Libraries for the solver
set(lib_vtk ${SV_LIB_THREEDSOLVER_VTK_NAME}${SV_MPI_NAME_EXT})
set(lib_leslib ${SV_LIB_THREEDSOLVER_DUMMY_LESLIB_NAME}${SV_MPI_NAME_EXT})
set(lib_restart ${SV_LIB_THREEDSOLVER_WRITE_RESTART_NAME}${SV_MPI_NAME_EXT})
set(lib_flowsolver ${SV_LIB_THREEDSOLVER_FORTRAN_NAME}${SV_MPI_NAME_EXT})
# Include files for the USE_VTK VERSION
if (SV_USE_VTK)
  set(CXXSRCS ${CXXSRCS} use_vtk.cxx)
  add_library(${lib_vtk} ${SV_LIBRARY_TYPE} use_vtk.cxx)
  target_link_libraries(${lib_vtk}
    ${VTK_LIBRARIES})
endif()

if(SV_USE_LESLIB)
  include_directories(${LESLIB_INCLUDE_DIR})
else()# dummy leslib
  add_library(${lib_leslib} ${SV_LIBRARY_TYPE} leslib_dummy_prototypes.c)
  set(LESLIB_LIBRARY ${lib_leslib})
endif()
add_library(${lib_restart} ${SV_LIBRARY_TYPE} ${CSRCS})
target_link_libraries(${lib_restart} ${SV_LIB_SOLVERIO_NAME}
  ${lib_leslib})
if(SV_USE_DUMMY_MPI)
  target_link_libraries(${lib_restart}
    ${SV_LIB_DUMMY_MPI_NAME} ${SV_LIB_DUMMY_FORTRAN_MPI_NAME})
else()
  target_link_libraries(${lib_restart}
    ${SV_MPI_EXTRA_LIBRARY} ${SV_MPI_Fortran_LIBRARIES})
endif()

add_library(${lib_flowsolver} ${SV_LIBRARY_TYPE} ${FSRCS})
#set_TARGET_PROPERTIES(${lib_flowsolver} PROPERTIES COMPILE_FLAGS "-extend_source 132")
  target_link_libraries(${lib_flowsolver}
    ${SV_LIB_SOLVERIO_NAME} ${lib_restart})
if(SV_USE_SVLS)
  target_link_libraries(${lib_flowsolver} ${SV_LIB_SVLS_NAME}${SV_MPI_NAME_EXT})
endif()
if(SV_USE_VTK)
  target_link_libraries(${lib_flowsolver} ${lib_vtk})
endif()
if(SV_USE_DUMMY_MPI)
  target_link_libraries(${lib_flowsolver} ${SV_LIB_DUMMY_MPI_NAME}
    ${SV_LIB_DUMMY_FORTRAN_MPI_NAME})
else()
  target_link_libraries(${lib_flowsolver} ${SV_MPI_EXTRA_LIBRARY}
    ${SV_MPI_Fortran_LIBRARIES})
endif()
set(HDRS Input.h ValType.h common_c.h)

simvascular_add_executable(${SV_FLOWSOLVER_EXE}${SV_MPI_NAME_EXT}
  SRCS ${CXXSRCS} ${CSRCS}
  DEV_SCRIPT_NAME "mysolver${SV_MPI_NAME_EXT}"
  INSTALL_SCRIPT_NAME "svsolver${SV_MPI_NAME_EXT}"
  INSTALL_COMP CoreExecutables
  INSTALL_DESTINATION ${SV_INSTALL_RUNTIME_DIR})

if(SV_USE_VTK)
  set(VTKLIB ${VTK_LIBRARIES})
else()
  set(VTKLIB "")
endif()

target_link_libraries(${SV_FLOWSOLVER_EXE}${SV_MPI_NAME_EXT} ${lib_flowsolver}
  ${GLOBAL_LIBRARIES} ${INTELRUNTIME_LIBRARIES}
  ${METIS_LIBRARY} ${SV_LIB_SVLS_NAME}${SV_MPI_NAME_EXT}
  ${LESLIB_LIBRARY}
  ${NSPCG_LIBRARY} ${SV_LIB_SOLVERIO_NAME}
  ${VTKLIB}
  )
if(SV_USE_DUMMY_MPI)
  target_link_libraries(${SV_FLOWSOLVER_EXE}${SV_MPI_NAME_EXT} ${SV_LIB_DUMMY_MPI_NAME}
    ${SV_LIB_DUMMY_FORTRAN_MPI_NAME})
else()
  target_link_libraries(${SV_FLOWSOLVER_EXE}${SV_MPI_NAME_EXT} ${SV_MPI_EXTRA_LIBRARY})
endif()

if(SV_INSTALL_LIBS)
  install(TARGETS ${lib_flowsolver} ${lib_restart} ${lib_leslib}
    RUNTIME DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT CoreExecutables
    LIBRARY DESTINATION ${SV_INSTALL_LIBRARY_DIR} COMPONENT CoreLibraries
    ARCHIVE DESTINATION ${SV_INSTALL_ARCHIVE_DIR} COMPONENT CoreLibraries
    )
  if(SV_USE_VTK)
    install(TARGETS ${lib_vtk}
      RUNTIME DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT CoreExecutables
      LIBRARY DESTINATION ${SV_INSTALL_LIBRARY_DIR} COMPONENT CoreLibraries
      ARCHIVE DESTINATION ${SV_INSTALL_ARCHIVE_DIR} COMPONENT CoreLibraries
      )
  endif()
endif()
if(SV_INSTALL_HEADERS)
  install(FILES ${HDRS}
    DESTINATION ${SV_INSTALL_INCLUDE_DIR}/core COMPONENT CoreHeaders
    )
endif()



