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

set(lib ${SV_LIB_SVFSILS_NAME}${SV_MPI_NAME_EXT})

set(FSRCS ADDBCMUL.f
  BCAST.f
  CGRAD.f
  CPUT.f
  FSISOLVER.f
  GMRES.f
  LHS.f
  NORM.f
  OMPLA.f
  PRECOND.f
  SPARMUL.f
  BC.f
  BICGS.f
  COMMU.f
  DOT.f
  GE.f
  INCOMMU.f
  LS.f
  NSSOLVER.f
  PCGMRES.f
  SOLVE.f)

set(CSRCS HRCPUT.c)

add_library(${lib} ${SV_LIBRARY_TYPE} ${FSRCS} ${CSRCS})

target_link_libraries(${lib} ${MPI_LIBRARY} ${MPI_Fortran_LIBRARIES})

# extra MPI libraries only if there are not set to NOT_FOUND or other null
if(SV_MPI_EXTRA_LIBRARY)
  target_link_libraries(${lib} ${SV_MPI_EXTRA_LIBRARY})
endif()

if(SV_INSTALL_LIBS)
  install(TARGETS ${lib}
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
