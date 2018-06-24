include(${CMAKE_CURRENT_SOURCE_DIR}/SimVascularScriptMacros.cmake)
file(MAKE_DIRECTORY ${TEMP_DIR}/developer ${TEMP_DIR}/install)

set(core_str)
if(UNIX)
  file(READ ${SV_SOURCE_DIR}/Scripts/home-locate-linux.sh locate_script)
  set(core_str "${locate_script}\n")
  set_env_string_concat(SV_HOME "$DIR" core_str)
endif()
if(WIN32 AND NOT UNIX)
  set(core_str "@echo off\r\n")
  set_env_string_concat(SV_HOME "%~dp0" core_str)
endif()

set(core_vars 	SV_OS SOFTWARE_PLATFORM COMPILER_VERSION
  SV_VERSION OPEN_SOFTWARE_BINARIES_TOPLEVEL
  LICENSED_SOFTWARE_TOPLEVEL SV_STATIC_BUILD)


foreach(var ${core_vars})
  set_env_string_concat(${var} "\@${var}\@" core_str)
  string(CONFIGURE "${core_str}" core_str @ONLY)
endforeach()


set(SCRIPT_FILES ${HOME_FILE_STRING})
if(UNIX)
  set(WIN_PATH )
  set(WIN_EXE )
  set(WIN_BAT )
  set(OSTYPE linux)
endif()

#Generate output scripts
if(WIN32 AND NOT CYGWIN)
  set(WIN_EXE ".exe")
  set(WIN_BAT ".bat")
  set(OSTYPE windows)
endif()

set(developer_script_string)
set(install_script_string)
set(FRAGMENT_INJECTION_STRING)
#---
# Write extra FILES
env_variable_to_value_variable(ENV_SV_HOME SV_HOME)

#configure_file(${SV_SOURCE_DIR}/../Distribution/windows/registry-patch-notargets.xml.in "${TEMP_DIR}/patch.xml" @ONLY)
set(PREAMBLE_DEVELOPER_FILE "${TEMP_DIR}/developer-preample")
set(PREAMBLE_INSTALL_FILE "${TEMP_DIR}/install-preample")
set(EXE_FILE_STRING "")

#----------------------------------------------------------------------------
# add runtime and lib paths
if(WIN32)
  set(developer_script_string "${developer_script_string}\nREM svSolver runtime and library locations\n")
  set(install_script_string "${install_script_string}\nREM svSolver runtime and library locations\n")
else()
  set(developer_script_string "${developer_script_string}\n# svSolver runtime and library locations\n")
  set(install_script_string "${install_script_string}\n# svSolver runtime and library locations\n")
endif()
append_env_string_concat(${ENV_LIBRARY_PATH_VARIABLE} "${ENV_SV_HOME}/${SV_INSTALL_LIBRARY_DIR}" install_script_string)
append_env_string_concat(${ENV_LIBRARY_PATH_VARIABLE} "${ENV_SV_HOME}/${SV_INSTALL_HOME_DIR}" install_script_string)
append_env_string_concat(${ENV_LIBRARY_PATH_VARIABLE} "${ENV_SV_HOME}/${SV_INSTALL_RUNTIME_DIR}" install_script_string)

append_env_string_concat(${ENV_LIBRARY_PATH_VARIABLE} "${CMAKE_BINARY_DIR}/bin" developer_script_string)
append_env_string_concat(${ENV_LIBRARY_PATH_VARIABLE} "${CMAKE_BINARY_DIR}/lib" developer_script_string)
append_env_string_concat(${ENV_LIBRARY_PATH_VARIABLE} "${CMAKE_BINARY_DIR}/bin/RelWithDebInfo" developer_script_string)
#----------------------------------------------------------------------------


file(WRITE ${PREAMBLE_DEVELOPER_FILE} "${core_str}\n")
file(APPEND ${PREAMBLE_DEVELOPER_FILE} "${developer_script_string}\n")

file(WRITE ${PREAMBLE_INSTALL_FILE} "${core_str}\n")
file(APPEND ${PREAMBLE_INSTALL_FILE} "${install_script_string}\n")

set(SCRIPT_FILES "${TEMP_DIR}/\@mode\@-preample")

# CHANGE FOR EXECUTABLE RENAME REMOVE
set(SCRIPT_FILES ${SCRIPT_FILES} "${SV_SOURCE_DIR}/Scripts/executable/\@mode\@/\@tar\@-\@OSTYPE\@.in")

unset(GDB)
unset(GDB_FLAGS)
set(INSTALL_SYMLINK_FILE "${TEMP_DIR}/setup-symlinks.sh")
file(WRITE ${INSTALL_SYMLINK_FILE} "${locate_script}\n\necho $DIR\n\n")
dev_message("SV_SCRIPT_TARGETS: ${SV_SCRIPT_TARGETS}")
foreach(tar ${SV_SCRIPT_TARGETS})
  simvascular_configure_script(
    TARGET "${tar}"
    FILES ${SCRIPT_FILES}
    INSTALL
    )
  if(WIN32 AND NOT CYGWIN)
    message("${tar}")
    ADD_CUSTOM_TARGET(${tar}-copy ALL
      DEPENDS ${tar}
      COMMAND ${CMAKE_COMMAND} -E remove ${OUTBIN_DIR}/$<TARGET_FILE_NAME:${tar}>
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${tar}> ${OUTBIN_DIR}
      COMMENT "Copying ${tar} to ${OUTBIN_DIR}"
      )
  endif()
  if(UNIX)
    file(APPEND ${INSTALL_SYMLINK_FILE} "rm -f /usr/local/bin/${${tar}_INSTALL_SCRIPT_NAME}${WIN_BAT}\n")
    file(APPEND ${INSTALL_SYMLINK_FILE} "ln -s $DIR/${${tar}_INSTALL_SCRIPT_NAME}${WIN_BAT} /usr/local/bin/${${tar}_INSTALL_SCRIPT_NAME}${WIN_BAT}\n")
    #file(APPEND ${INSTALL_SYMLINK_FILE} "chmod a+x /usr/local/bin/${${tar}_INSTALL_SCRIPT_NAME}${WIN_BAT}\n")
  endif()
endforeach()
if(UNIX)
  install(PROGRAMS ${INSTALL_SYMLINK_FILE} DESTINATION ${SV_INSTALL_SCRIPT_DIR})
  install(FILES "README-install" DESTINATION ${SV_INSTALL_SCRIPT_DIR})
  #install(CODE "execute_process(COMMAND \"${INSTALL_SYMLINK_FILE}\"")")
endif()
if(BUILD_TESTING)
  configure_file(executable_names.tcl.in ${SV_TEST_DIR}/common/executable_names.tcl)
endif()
