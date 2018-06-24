
set(lib ${SV_LIB_SVLS_NAME}${SV_MPI_NAME_EXT})

set(FSRCS ADDBCMUL.f
  BC.f
  #BICGS.f
  CGRAD.f
  COMMU.f
  CPUT.f
  DOT.f
  GE.f
  GMRES.f
  INCOMMU.f
  LHS.f
  LS.f
  NORM.f
  NSSOLVER.f
  #OMPLA.f
  PRECOND.f
  SOLVE.f
  SPARMUL.f)


add_library(${lib} ${SV_LIBRARY_TYPE} ${FSRCS})

if(SV_USE_DUMMY_MPI)
  target_link_libraries(${lib} ${SV_LIB_DUMMY_FORTRAN_MPI_NAME})
else()
  target_link_libraries(${lib} ${MPI_LIBRARY} ${MPI_Fortran_LIBRARIES})
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
