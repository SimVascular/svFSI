include_directories(./)

set(CSRCS comm.c util.c debug.c setup.c grsetup.c timer.c
       node_refine.c initmsection.c
       xyzpart.c pspases.c frename.c
       iintsort.c iidxsort.c ikeysort.c ikeyvalsort.c
       kmetis.c gkmetis.c ometis.c
       initpart.c match.c
       kwayfm.c kwayrefine.c kwaybalance.c
       remap.c stat.c fpqueue.c
       ametis.c rmetis.c initbalance.c
       mdiffusion.c diffutil.c wave.c
       csrmatch.c redomylink.c balancemylink.c
       selectq.c akwayfm.c serial.c move.c
       mmetis.c mesh.c memory.c weird.c backcompat.c)

if(SV_USE_THIRDPARTY_SHARED_LIBRARIES)
  add_library(${PARMETIS_SVFSI_LIBRARY_NAME} SHARED ${CSRCS})
else()
  add_library(${PARMETIS_SVFSI_LIBRARY_NAME} STATIC ${CSRCS})
endif()

target_link_libraries(${PARMETIS_SVFSI_LIBRARY_NAME} ${METIS_SVFSI_LIBRARY_NAME})

if(SV_INSTALL_LIBS)
  install(TARGETS ${PARMETIS_SVFSI_LIBRARY_NAME}
  RUNTIME DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT ThirdPartyExecutables
  LIBRARY DESTINATION ${SV_INSTALL_LIBRARY_DIR} COMPONENT ThirdPartyLibraries
  ARCHIVE DESTINATION ${SV_INSTALL_ARCHIVE_DIR} COMPONENT ThirdPartyLibraries)
endif()

if(SV_INSTALL_HEADERS)

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/temp)
  set(TEMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/temp)
  set(HDRS)
  #copy and rename files
  file(READ "parmetislib.h" METIS_H_TEMP)

  set(filename "macros.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "defs.h")
  file(COPY ${filename} DESTINATION ${TEMP_DIR})
  file(RENAME ${TEMP_DIR}/${filename}
       ${TEMP_DIR}/metis_${filename} )
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "proto.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "rename.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "stdheaders.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "struct.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "parmetis.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename ${TEMP_DIR}/parmetislib.h)
  file(WRITE ${filename} ${METIS_H_TEMP})
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi
    COMPONENT ThirdPartyHeaders)

endif()
