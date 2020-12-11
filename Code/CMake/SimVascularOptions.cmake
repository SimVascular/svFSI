# Initial SV Options
#-----------------------------------------------------------------------------
# Developer flag (Output extra info during configure)
option(SV_DEVELOPER_OUTPUT "This is a developer mode to print extra messages during configure" OFF)

set(SV_BUILD_TYPE "CMAKE" CACHE STRING "Designate CMAKE build" FORCE)
set_property(CACHE SV_BUILD_TYPE PROPERTY STRINGS CMAKE)
mark_as_advanced(SV_BUILD_TYPE)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Distribution
option(SV_ENABLE_DISTRIBUTION "Distribute" OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Enable Testing
option(BUILD_TESTING "Build ${PROJECT_NAME} testing" OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Libs - SHARED or STATIC
option(BUILD_SHARED_LIBS "Build ${PROJECT_NAME} as shared libraries." OFF)

set(SV_LIBRARY_TYPE "STATIC" CACHE STRING "Options are STATIC or SHARED" FORCE)
set_property(CACHE SV_LIBRARY_TYPE PROPERTY STRINGS STATIC SHARED)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# SimVascular Build options
option(SV_SUPPRESS_WARNINGS "Option to suppress all compiler warnings while compiling" ON)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# General Options
option(SV_USE_MPI "Use MSMPI" ON)
option(SV_USE_MSMPI "Use MSMPI" OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# ThirdParty
option(SV_USE_ZLIB "Use ZLib" ON)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Remaining optional dependencies
#-----------------------------------------------------------------------------
# Enable Intel Runtime libs if we need or want them
option(SV_USE_INTEL "Add Intel Runtime Libraries (these may be needed by some libraries)" OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# All OS
option(SV_USE_NOTIMER "Use notimer" ON)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Solver Build Options (Modules)
option(SV_USE_METIS_SVFSI "Use metis_svfsi Library" ON)
option(SV_USE_PARMETIS_SVFSI "Use parmetis_svfsi Library" ON)
option(SV_USE_TETGEN "Use tetgen Library" ON)
option(SV_USE_TRILINOS "Use Trilinos Library with svFSI" OFF)
option(SV_USE_HYPRE "Use Hypre Library with svFSI" ON)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Externals
set(SV_EXTERNALS_TOPLEVEL_DIR "${CMAKE_BINARY_DIR}/sv_externals" CACHE PATH "Externals toplevel directory")
set(SV_EXTERNALS_SRC_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/src" CACHE PATH "Externals toplevel src dir")
set(SV_EXTERNALS_BLD_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/build" CACHE PATH "Externals toplevel build dir")
set(SV_EXTERNALS_PFX_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/prefix" CACHE PATH "Externals toplevel prefix dir")
set(SV_EXTERNALS_BIN_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/bin/${SV_COMPILER_DIR}/${SV_COMPILER_VERSION_DIR}/${SV_ARCH_DIR}" CACHE PATH "Externals toplevel bin dir")

set(SV_EXTERNALS_INSTALL_PREFIX "sv_externals" CACHE PATH "Externals toplevel directory")

option(SV_EXTERNALS_USE_TOPLEVEL_DIR "If ON, SV_EXTERNALS_TOPLEVEL_DIR will be used as location for external packages" OFF)


#-----------------------------------------------------------------------------
# SVFSILS linear solver is always on 
#-----------------------------------------------------------------------------
set(USE_SVFSILS 1)
set(SVFSILS_BUILD_TYPE "Source")

#-----------------------------------------------------------------------------
# WIN32
option(SV_USE_WIN32_REGISTRY "Use Windows registry to obtain certain settings (install mode)" OFF)
mark_as_advanced(SV_USE_WIN32_REGISTRY)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Enable Fortran
enable_language(Fortran)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # svsolver requires -ffixed-line-length-132 but
  # svFSI does not compile with this flag
  # reset flags for svFSILS and svFSI in their local CMakeLists
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132 -cpp")
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132 -fpp")
endif()
#-----------------------------------------------------------------------------
