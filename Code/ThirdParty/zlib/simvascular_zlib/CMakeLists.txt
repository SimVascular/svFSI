project(SimVascularZLIB)

include_directories(BEFORE ${SimVascularZLIB_SOURCE_DIR})

# source files for zlib
set(SRCS
  adler32.c  compress.c  crc32.c  deflate.c  gzio.c   inffast.c
  inflate.c  inftrees.c  trees.c  uncompr.c  zutil.c
  )
set(HDRS crc32.h deflate.h inffast.h inffixed.h inflate.h inftrees.h
  trees.h zconf.h zlib.h zutil.h)

if(WIN32)
  if(BUILD_SHARED_LIBS)
    set(ZLIB_DLL 1)
  endif(BUILD_SHARED_LIBS)
endif(WIN32)

set(lib ${SV_LIB_THIRDPARTY_ZLIB_NAME})
if(SV_USE_THIRDPARTY_SHARED_LIBRARIES)
  add_library(${lib} SHARED ${SRCS})
else()
  add_library(${lib} STATIC ${SRCS})
endif()

if(SV_INSTALL_LIBS)
  install(TARGETS ${lib}
    RUNTIME DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT ThirdPartyExecutables
    LIBRARY DESTINATION ${SV_INSTALL_LIBRARY_DIR} COMPONENT ThirdPartyLibraries
    ARCHIVE DESTINATION ${SV_INSTALL_ARCHIVE_DIR} COMPONENT ThirdPartyLibraries
    )
endif()
if(SV_INSTALL_HEADERS)
  install(FILES ${HDRS}
    DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/zlib COMPONENT ThirdPartyHeaders
    )
endif()
