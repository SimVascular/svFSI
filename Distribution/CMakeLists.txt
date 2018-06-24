#-----------------------------------------------------------------------------
# Include Macros
#-----------------------------------------------------------------------------
include(InstallRequiredSystemLibraries)
include(CPackComponent)

#-----------------------------------------------------------------------------
# Install program
#-----------------------------------------------------------------------------
install(PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
  DESTINATION ${SV_INSTALL_RUNTIME_DIR}
  COMPONENT CoreExecutables)

if(LINUX)
  set(SV_INSTALLER_TYPE DEB CACHE STRING "When distributing, specify the installer package type.")
    set_property(CACHE SV_INSTALLER_TYPE PROPERTY STRINGS DEB TGZ)
endif()

#-----------------------------------------------------------------------------
# Platform independent options
#-----------------------------------------------------------------------------
set(CPACK_PACKAGE_NAME "svSolver")
set(CPACK_PACKAGE_VENDOR "SimVascular.org")
set(CPACK_PACKAGE_CONTACT "SimVascular<sv@simvascular.org>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SimVascular Solvers")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license.txt")
set(CPACK_PACKAGE_VERSION_MAJOR ${SV_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${SV_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${SV_PATCH_VERSION})
set(INSTALL_TYPE_STR)
string(REPLACE "." "-" SV_INSTALL_VERSION_DIR ${SV_FULL_VERSION})

#-----------------------------------------------------------------------------
# WINDOWS
#-----------------------------------------------------------------------------
if(WIN32)
  # Package name
   set(CPACK_PACKAGE_FILE_NAME "svsolver-windows-win-x64-${SV_FULL_VERSION}${INSTALL_TYPE_STR}")

  # MPI dlls
  find_path(INTEL_REDIST_DIR names fmpich2.dll PATHS $ENV{PATH})
  file(GLOB INTEL_REDIST_DLLS "${INTEL_REDIST_DIR}/*.dll")
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS} ${INTEL_REDIST_DLLS})
  set(RM_DLL "fmpich2.dll;fmpich2g.dll;mpich2mpi.dll")
  set(RM_DLL_FULL)
  foreach(dll ${RM_DLL})
  	sv_list_match(tmp ${dll} ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS})
  	list(APPEND RM_DLL_FULL ${tmp})
  endforeach()
  list(REMOVE_ITEM CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${RM_DLL_FULL})

  # Cpack install help
  set(CPACK_MONOLITHIC_INSTALL 1)
  set(CPACK_COMPONENTS_IGNORE_GROUPS 1)
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS} ${INTEL_REDIST_DLLS})

  # Visual studios
  if(MSVC)
    set(CPACK_GENERATOR "WIX")
    set(CPACK_WIX_EXTENSIONS WixUtilExtension)
    set(CPACK_WIX_UI_BANNER ${CMAKE_CURRENT_SOURCE_DIR}/windows/msi-banner.jpg)
    set(CPACK_WIX_UI_DIALOG ${CMAKE_CURRENT_SOURCE_DIR}/windows/msi-dialog.jpg)
    if(SV_RELEASE_TYPE MATCHES "Release")
      set(CPACK_PACKAGE_INSTALL_DIRECTORY "SimVascular/Release")
      set(CPACK_WIX_PROGRAM_MENU_FOLDER "SimVascular")
      set(CPACK_WIX_PRODUCT_GUID "39068F8A-D4E9-405F-B22B-9B417849AEE6")
      set(CPACK_WIX_UPGRADE_GUID "2359EF69-C03F-401C-9BC3-C19B10297E3B")
      set(PACKAGE_NAME "svsolver")
    elseif(SV_RELEASE_TYPE MATCHES "Beta")
      set(CPACK_WIX_PROGRAM_MENU_FOLDER "SimVascular-Beta")
      set(PACKAGE_NAME "svsolver-beta")
      set(CPACK_PACKAGE_INSTALL_DIRECTORY "SimVascular/${TIMESTAMP}")
      set(CPACK_WIX_PRODUCT_GUID "9453BB7E-CD7E-46ED-AA9E-5514B56FC901")
      set(CPACK_WIX_UPGRADE_GUID "39C845FF-68E1-4AA7-B329-3407B47B3361")
    endif()
    set(CPACK_WIX_PATCH_FILE "${TEMP_DIR}/patch.xml")
    set(CPACK_WIX_PRODUCT_ICON "${CMAKE_CURRENT_SOURCE_DIR}/windows/SimVascular.ico")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/windows/mywix.template-notargets.in" "${CMAKE_CURRENT_BINARY_DIR}/mywix.template")
    set(CPACK_WIX_TEMPLATE "${CMAKE_CURRENT_BINARY_DIR}/mywix.template")
  endif()
endif()

#-----------------------------------------------------------------------------
# LINUX
#-----------------------------------------------------------------------------
if(LINUX)
  # Package name
  set(CPACK_PACKAGE_FILE_NAME "svsolver-linux-x64-${SV_FULL_VERSION}${INSTALL_TYPE_STR}")

    if(SV_INSTALLER_TYPE STREQUAL DEB)
        set(CPACK_GENERATOR "DEB")
        set(CPACK_PACKAGING_INSTALL_PREFIX "/usr/local/sv/svsolver/${SV_INSTALL_VERSION_DIR}")
    elseif(SV_INSTALLER_TYPE STREQUAL TGZ)
      set(CPACK_GENERATOR "TGZ")
      set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY ON)
    endif()
endif()

#-----------------------------------------------------------------------------
# APPLE
#-----------------------------------------------------------------------------
if(APPLE)
  # Package name
  set(CPACK_PACKAGE_FILE_NAME "svsolver-darwin-x64-${SV_FULL_VERSION}${INSTALL_TYPE_STR}")

  # Components
  set(CPACK_COMPONENTS_ALL CoreExecutables CoreLibraries CoreHeaders
    ThirdPartyLibraries ThirdPartyHeaders MPIExecutables MPILibraries
    ExecutableScripts)

  # Cpack install help
  set(CPACK_GENERATOR "productbuild")
  #configure_file(${CMAKE_CURRENT_SOURCE_DIR}/OSX/Info.plist.in
  #  ${CMAKE_CURRENT_BINARY_DIR}/OSX/Info.plist @ONLY)

  # Install location
  set(CPACK_PACKAGING_INSTALL_PREFIX "/usr/local/sv/svsolver/${SV_INSTALL_VERSION_DIR}")
endif()

if(CPACK_PACKAGE_FILENAME_TAG)
	set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-${CPACK_PACKAGE_FILENAME_TAG}")
endif()

message(STATUS "Packname name: ${CPACK_PACKAGE_FILE_NAME}")
include(CPack)
