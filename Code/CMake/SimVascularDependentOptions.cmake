# Options and defined depending on the system, and the intial options
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# APPLE
if(APPLE)
  set(CMAKE_OSX_ARCHITECTURES "" CACHE STRING "" FORCE)
  # Note: By setting CMAKE_OSX_* variables before any enable_language() or project() calls,
  #       we ensure that the bitness will be properly detected.
  set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/opt/local/lib")
endif()

#-----------------------------------------------------------------------------
# WIN32
if(WIN32)
  set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -DWINDOWS -DWIN32")
  if(NOT IS64)
    if(NOT "${CMAKE_GENERATOR}" MATCHES ".*Win64")
      set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -D_X86_")
    endif()
  endif()
  if(CYGWIN)
    set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -DCYGWIN")
  endif()
  set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -DSV_WRAP_FORTRAN_IN_CAPS_NO_UNDERSCORE")
  set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -D__VC__")
  check_library_exists("${CMAKE_CXX_STANDARD_LIBRARIES}" gethostname "" HAVE_STDGETHOSTNAME)
  if(NOT HAVE_STDGETHOSTNAME)
    check_library_exists("wsock32.lib" gethostname "" HAVE_WSOCK_GETHOSTNAME)
    if(HAVE_WSOCK_GETHOSTNAME)
      set (CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} wsock32.lib")
    else()
      MESSAGE(AUTHOR_WARNING "gethostname has not beed found! The flowsolver will not compile")
    endif()
  endif()

  set (CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} Shlwapi.lib")
endif()

#-----------------------------------------------------------------------------
# LINUX
if(UNIX)
  set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -DUNIX")
  set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -DSV_WRAP_FORTRAN_IN_LOWERCASE_WITH_UNDERSCORE")
endif()
if(LINUX)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pthread -static")
endif()

#-----------------------------------------------------------------------------
# Visual Studio flags
if(MSVC)
  set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -DMSVC /EHsc")
# SUPPRESS_VC_DEPRECATED_WARNINGS
  set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS")
else()
  if (BUILD_SHARED_LIBS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
endif()

# Compiler Flags
#-----------------------------------------------------------------------------
if(SV_SUPPRESS_WARNINGS)
  add_definitions("-w")
endif()
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpermissive")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
  set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -DGCC")
endif()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pthread")

#-----------------------------------------------------------------------------
# Set a default build type (if none was specified)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  mark_as_advanced(CMAKE_BUILD_TYPE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
  	"MinSizeRel" "RelWithDebInfo")
endif()

#-----------------------------------------------------------------------------
# Shared Libs
#-----------------------------------------------------------------------------
set(SV_INSTALL_HEADERS ON)
set(SV_INSTALL_EXTERNALS ON)
set(SV_INSTALL_LIBS ON)
if(BUILD_SHARED_LIBS)
  set(SV_LIBRARY_TYPE "SHARED" CACHE STRING "Shared cache" FORCE)
  set(SV_STATIC_BUILD "0")
else()
  set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -DSV_STATIC_LINK -DSV_STATIC_BUILD")
  set(SV_STATIC_BUILD "1")
endif()

#-----------------------------------------------------------------------------
# Flowsolver
#-----------------------------------------------------------------------------
if(SV_SOLVERIO_REDIRECT)
  set(GLOBAL_DEFINES "${GLOBAL_DEFINES} -DBUILD_WITH_FLOWSOLVER_STDOUT_STDERR_REDIRECT")
endif()

#-----------------------------------------------------------------------------
# Postsolver
#-----------------------------------------------------------------------------
if (NOT SV_USE_VTK)
  set(SV_USE_SVPOST "OFF" CACHE BOOL "Cannot build svpost without vtk" FORCE)
  set(SV_USE_SVPRE "OFF" CACHE BOOL "Cannot build svpre without vtk" FORCE)
endif()
