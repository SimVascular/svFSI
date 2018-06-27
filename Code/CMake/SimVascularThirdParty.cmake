#-----------------------------------------------------------------------------
# ZLIB
if(SV_USE_ZLIB)
  set(USE_ZLIB ON)
  simvascular_third_party(zlib)
  if(NOT SV_USE_SYSTEM_ZLIB)
    set(ZLIB_LIBRARY ${SV_LIB_THIRDPARTY_ZLIB_NAME})
  else()
    find_package(ZLIB)
  endif()
else()
  unset(ZLIB_LIBRARY CACHE)
endif()

#-----------------------------------------------------------------------------
# METIS_SVFSI
if(SV_USE_METIS_SVFSI)
  set(USE_METIS_SVFSI ON)
  simvascular_third_party(metis_svfsi)
  # require to be built here
  set(METIS_SVFSI_LIBRARY ${SV_LIB_THIRDPARTY_METIS_SVFSI_NAME})
endif()

#-----------------------------------------------------------------------------
# PARMETIS_SVFSI
if(SV_USE_PARMETIS_SVFSI)
  set(USE_PARMETIS_SVFSI ON)
  simvascular_third_party(parmetis_svfsi)
  # require to be built here 
  set(PARMETIS_SVFSI_LIBRARY ${SV_LIB_THIRDPARTY_PARMETIS_SVFSI_NAME})
endif()

#-----------------------------------------------------------------------------
# TETGEN
if(SV_USE_TETGEN)
  set(USE_TETGEN ON)
  simvascular_third_party(tetgen)
  # require to be built here 
  set(TETGEN_LIBRARY ${SV_LIB_THIRDPARTY_TETGEN_NAME})
endif()
