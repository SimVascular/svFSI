
# reset flags if using gnu compilers since -ffixed-line-length-132 does not work here 
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-cpp -pthread -DWITH_TRILINOS")
else()
  # nothing for now
  # may need to reset for intel compiler or others
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
