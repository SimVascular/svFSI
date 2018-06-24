include_directories(./)

set(CSRCS coarsen.c fm.c initpart.c match.c ccgraph.c 
       pmetis.c pqueue.c refine.c util.c timing.c debug.c 
       bucketsort.c graph.c stat.c kmetis.c kwayrefine.c 
       kwayfm.c balance.c ometis.c srefine.c sfm.c separator.c 
       mincover.c mmd.c mesh.c meshpart.c frename.c fortran.c 
       myqsort.c compress.c parmetis.c estmem.c 
       mpmetis.c mcoarsen.c mmatch.c minitpart.c mbalance.c 
       mutil.c mkmetis.c mkwayrefine.c mkwayfmh.c 
       mrefine2.c minitpart2.c mbalance2.c mfm2.c 
       kvmetis.c kwayvolrefine.c kwayvolfm.c subdomains.c 
       mfm.c memory.c mrefine.c checkgraph.c)

if(SV_USE_THIRDPARTY_SHARED_LIBRARIES)
  add_library(${METIS_SVFSI_LIBRARY_NAME} SHARED ${CSRCS})
else()
  add_library(${METIS_SVFSI_LIBRARY_NAME} STATIC ${CSRCS})
endif()

if(SV_INSTALL_LIBS)
  install(TARGETS ${METIS_SVFSI_LIBRARY_NAME}
  RUNTIME DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT ThirdPartyExecutables
  LIBRARY DESTINATION ${SV_INSTALL_LIBRARY_DIR} COMPONENT ThirdPartyLibraries
  ARCHIVE DESTINATION ${SV_INSTALL_ARCHIVE_DIR} COMPONENT ThirdPartyLibraries)
endif()

if(SV_INSTALL_HEADERS)

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/temp)
  set(TEMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/temp)
  set(HDRS)
  #copy and rename files
  file(READ "metis.h" METIS_H_TEMP)

  set(filename "macros.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "defs.h")
  file(COPY ${filename} DESTINATION ${TEMP_DIR})
  file(RENAME ${TEMP_DIR}/${filename}
       ${TEMP_DIR}/metis_${filename} )
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "proto.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "rename.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "stdheaders.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "struct.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "parmetis.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename ${TEMP_DIR}/metis.h)
  file(WRITE ${filename} ${METIS_H_TEMP})
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi
    COMPONENT ThirdPartyHeaders)
endif()


