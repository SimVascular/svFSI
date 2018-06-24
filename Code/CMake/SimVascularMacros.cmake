# Copyright (c) 2014-2015 The Regents of the University of California.
# All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

include (CMakeParseArguments)

macro(dev_message string)
  if(SV_DEVELOPER_OUTPUT)
    message("DEV: ${string}")
  endif()
endmacro()

#-----------------------------------------------------------------------------
# simvascular_external - Macro to find libraries needed by simvascular
# and create the necessary variables to load and link them.
#
macro(simvascular_external _pkg)
  string(TOLOWER "${_pkg}" _lower)

  dev_message("Configuring ${_pkg}")

  set(options OPTIONAL VERSION_EXACT
    DOWNLOADABLE SYSTEM_DEFAULT
    SVEXTERN_CONFIG ADD_INSTALL SHARED_LIB NO_MODULE
    )
  set(oneValueArgs VERSION)
  set(multiValueArgs PATHS HINTS COMPONENTS)

  cmake_parse_arguments("simvascular_external"
    "${options}"
    "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set(EXTRA_ARGS)
  if(simvascular_external_COMPONENTS)
    set(EXTRA_ARGS COMPONENTS ${simvascular_external_COMPONENTS})
  endif()

  if(simvascular_external_NO_MODULE)
    set(EXTRA_ARGS ${EXTRA_ARGS} NO_MODULE)
  endif()

  set(${_pkg}_VERSION ${simvascular_external_VERSION})
  if(simvascular_external_VERSION_EXACT)
    set(${_pkg}_VERSION ${${_pkg}_VERSION} EXACT)
  endif()

  unset(ARG_STRING)
  set(_paths "${simvascular_external_PATHS}")
  if(NOT simvascular_external_PATHS)
    set(_paths "${CMAKE_MODULE_PATH}")
  endif()

  #message(STATUS "Search paths for ${_pkg}Config.cmake: ${_paths}")

  if(simvascular_external_SYSTEM_DEFAULT)
    option(SV_USE_SYSTEM_${_pkg} "Use system ${_pkg}" ON)
    mark_as_advanced(SV_USE_SYSTEM_${_pkg})
  else()
    option(SV_USE_SYSTEM_${_pkg} "Use system ${_pkg}" OFF)
  endif()

  if((NOT SV_SUPERBUILD AND simvascular_external_SVEXTERN_CONFIG) OR
    (simvascular_external_SVEXTERN_CONFIG AND SV_USE_SYSTEM_${_pkg}))

  find_package(${_pkg} ${EXTRA_ARGS}
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/CMake
    NO_CMAKE_MODULE_PATH
    NO_DEFAULT_PATH)
  elseif(NOT SV_SUPERBUILD)
    find_package(${_pkg} ${EXTRA_ARGS})
  elseif(SV_USE_SYSTEM_${_pkg})
    find_package(${_pkg} ${EXTRA_ARGS})
  endif()

  if(simvascular_external_DOWNLOADABLE)
    set(SV_DEPENDS ${SV_DEPENDS} ${_pkg})
    list( REMOVE_DUPLICATES SV_DEPENDS )
  endif()

  if(SV_USE_${_pkg})
    set(USE_${_pkg} ON)
  endif()

  if(simvascular_external_SHARED_LIB)
    set(SV_EXTERNAL_SHARED_LIBS ${SV_EXTERNAL_SHARED_LIBS} ${_pkg})
  endif()

  if(${_pkg}_FOUND)
    message(STATUS "PKG ${_pkg} found!")
    if( ${_pkg}_INCLUDE_DIR )
      dev_message("Including dir: ${${_pkg}_INCLUDE_DIR}")
      # This get many of them
      include_directories(${${_pkg}_INCLUDE_DIR})
    endif()
    if(SV_INSTALL_EXTERNALS)
      if(simvascular_external_ADD_INSTALL)
      	getListOfVars("${_pkg}" "LIBRARY" ${_pkg}_VARS_INSTALL)
      	foreach(lib_install ${${_pkg}_VARS_INSTALL})
      	  list(APPEND ${_pkg}_LIBRARY_INSTALL "${${lib_install}}")
      	endforeach()
      	#message(STATUS "${_pkg}_LIBRARY_INSTALL: ${${_pkg}_LIBRARY_INSTALL}")
      	#install(FILES "${${_pkg}_LIBRARY_INSTALL}" DESTINATION ${SV_INSTALL_EXTERNAL_LIBRARY_DIR})
      endif()
    endif()
  endif()
  unset(simvascular_external_SVEXTERN_CONFIG)
  unset(simvascular_external_ADD_INSTALL)
  if(SV_DEVELOPER_OUTPUT)
    message(STATUS "Finished Configuring ${_pkg}")
    message(STATUS "")
  endif()
endmacro()
#-----------------------------------------------------------------------------
# unset_simvascular_external
#
macro(unset_simvascular_external _pkg)
  string(TOLOWER "${_pkg}" _lower)

  set(options OPTIONAL VERSION_EXACT DOWNLOADABLE SVEXTERN_DEFAULT SVEXTERN_CONFIG)
  set(oneValueArgs VERSION)
  set(multiValueArgs PATHS HINTS)

  cmake_parse_arguments("simvascular_external"
    "${options}"
    "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  unset(SV_USE_SYSTEM_${_pkg})
  list(REMOVE_ITEM SV_DEPENDS ${_pkg})
endmacro()

#-----------------------------------------------------------------------------
# simvascular_third_party
#
macro(simvascular_third_party _pkg)
  string(TOLOWER "${_pkg}" _lower)
  string(TOUPPER "${_pkg}" _upper)

  set(options OPTIONAL VERSION_EXACT
    DOWNLOADABLE SYSTEM_DEFAULT
    SVEXTERN_CONFIG ADD_INSTALL
    )
  set(oneValueArgs VERSION)
  set(multiValueArgs PATHS HINTS COMPONENTS)

  cmake_parse_arguments("simvascular_third_party"
    "${options}"
    "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  set(${_upper}_SUBDIR ThirdParty/${_pkg})
  if(simvascular_third_party_SYSTEM_DEFAULT)
    option(SV_USE_SYSTEM_${_upper} "Use system ${_pkg}" ON)
  else()
    option(SV_USE_SYSTEM_${_upper} "Use system ${_pkg}" OFF)
  endif()

  mark_as_advanced(SV_USE_SYSTEM_${_upper})

  configure_file(${SV_SOURCE_DIR}/${${_upper}_SUBDIR}/simvascular_${_lower}.h.in
    ${SV_BINARY_DIR}/${${_upper}_SUBDIR}/simvascular_${_lower}.h)
  include_directories(BEFORE ${SV_BINARY_DIR}/${${_upper}_SUBDIR} ${SV_SOURCE_DIR}/${${_upper}_SUBDIR})
  if(SV_USE_SYSTEM_${_upper})
    set(${_upper}_LIBRARIES)
    set(${_upper}_LIBRARY)
  else()
    if(NOT SV_SUPERBUILD)
      set(${_upper}_LIBRARY_NAME _simvascular_thirdparty_${_lower})
      add_subdirectory(${${_upper}_SUBDIR}/simvascular_${_lower})
    endif()
  endif()
endmacro()
#-----------------------------------------------------------------------------
# print_vars - THis is a simple marco to print out a list of variables
# with their names and value, used mostly for debugging
#
macro(print_vars _VARLIST)
  foreach(var ${${_VARLIST}})
    message(STATUS "${var}: ${${var}}")
  endforeach()
  message(STATUS "")
endmacro()

macro(dev_print_vars _VARLIST)
  foreach(var ${${_VARLIST}})
    dev_message("${var}: ${${var}}")
  endforeach()
  message(STATUS "")
endmacro()


#-----------------------------------------------------------------------------
# listvars2vars - THis is a simple marco to print out a list of variables
# with their names and value, used mostly for debugging
#
macro(listvars2vals _VARLIST _VALS)
  foreach(var ${${_VARLIST}})
    set(_vals ${_vals} ${${var}})
  endforeach()
  list(REMOVE_DUPLICATES _vals)
  set(${_VALS} ${_vals})
  unset(_vals)
endmacro()

#-----------------------------------------------------------------------------
# check_library_exists_concat -
#
macro(check_library_exists_concat LIBRARY SYMBOL VARIABLE)
  check_library_exists("${LIBRARY};${LINK_LIBS}" ${SYMBOL} "" ${VARIABLE})
  if(${VARIABLE})
    set(LINK_LIBS ${LINK_LIBS} ${LIBRARY})
  endif(${VARIABLE})
endmacro ()

#-----------------------------------------------------------------------------
# simvascular_add_executable -
macro(simvascular_add_executable TARGET_NAME)
  set(options NO_SCRIPT)
  set(oneValueArgs DEV_SCRIPT_NAME INSTALL_SCRIPT_NAME INSTALL_COMP INSTALL_DESTINATION)
  set(multiValueArgs SRCS)

  unset(simvascular_add_executable_INSTALL_SCRIPT_NAME)
  unset(simvascular_add_executable_DEV_SCRIPT_NAME)
  unset(simvascular_add_executable_NO_SCRIPT)

  cmake_parse_arguments("simvascular_add_executable"
    "${options}"
    "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set(WINDOWS_ICON_RESOURCE_FILE "")
  if(WIN32)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/icons/${TARGET_NAME}.rc")
      set(WINDOWS_ICON_RESOURCE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/icons/${TARGET_NAME}.rc")
    endif()
  endif()

  set(_app_compile_flags )
  if(WIN32)
    set(_app_compile_flags "${_app_compile_flags} -DPOCO_NO_UNWINDOWS -DWIN32_LEAN_AND_MEAN")
  endif()

  if(WIN32)
    add_executable(${TARGET_NAME} WIN32 ${simvascular_add_executable_SRCS} ${WINDOWS_ICON_RESOURCE_FILE})
  else()
    add_executable(${TARGET_NAME} ${simvascular_add_executable_SRCS} ${WINDOWS_ICON_RESOURCE_FILE})
  endif()

  set_target_properties(${TARGET_NAME} PROPERTIES
    COMPILE_FLAGS "${_app_compile_flags}")

  if(simvascular_add_executable_NO_SCRIPT)
    if(	simvascular_add_executable_DEV_SCRIPT_NAME OR simvascular_add_executable_INSTALL_SCRIPT_NAME )
      message(ERROR "Cannot specify no script and specify script names!")
    endif()
    set(${TARGET_NAME}_EXECUTABLE_NAME ${TARGET_NAME} CACHE INTERNAL "" FORCE)
  endif()
  if(NOT simvascular_add_executable_NO_SCRIPT)
    if(simvascular_add_executable_DEV_SCRIPT_NAME)
      set(SV_SCRIPT_TARGETS_WORK ${SV_SCRIPT_TARGETS})
      list(APPEND SV_SCRIPT_TARGETS_WORK "${TARGET_NAME}")
      list(REMOVE_DUPLICATES SV_SCRIPT_TARGETS_WORK)
      set(SV_SCRIPT_TARGETS ${SV_SCRIPT_TARGETS_WORK} CACHE INTERNAL "" FORCE)
      set(${TARGET_NAME}_DEVELOPER_SCRIPT_NAME ${simvascular_add_executable_DEV_SCRIPT_NAME} CACHE INTERNAL "" FORCE)
      set(${TARGET_NAME}_EXECUTABLE_NAME ${${TARGET_NAME}_DEVELOPER_SCRIPT_NAME} CACHE INTERNAL "" FORCE)
    endif()
    if(simvascular_add_executable_INSTALL_SCRIPT_NAME)
      set(${TARGET_NAME}_INSTALL_SCRIPT_NAME ${simvascular_add_executable_INSTALL_SCRIPT_NAME} CACHE INTERNAL "" FORCE)
    endif()
  endif()

  # CHANGE FOR EXECUTABLE RENAME REMOVE (re enable if statement)
  if(simvascular_add_executable_INSTALL_DESTINATION)
    if(simvascular_add_executable_INSTALL_COMP)
      set(_COMPARGS COMPONENT ${simvascular_add_executable_INSTALL_COMP})
    endif()
    if(APPLE)
      install(TARGETS ${TARGET_NAME}
        RUNTIME DESTINATION ${simvascular_add_executable_INSTALL_DESTINATION}
        ${_COMPARGS})
    else()
      install(TARGETS ${TARGET_NAME}
        RUNTIME DESTINATION ${simvascular_add_executable_INSTALL_DESTINATION}
        ${_COMPARGS})
    endif()
  endif()

endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
macro(simvascular_get_external_path_from_include_dir pkg)
  if(NOT ${pkg}_DIR)
    message(FATAL_ERROR "${pkg}_DIR is not set and must be set if using system")
  endif()
  if(${pkg}_INCLUDE_DIRS)
    list(GET ${pkg}_INCLUDE_DIRS 0 TMP_DIR)
    get_filename_component(TMP_DIR ${TMP_DIR} PATH)
    get_filename_component(TMP_DIR ${TMP_DIR} PATH)
  endif()
  if(${pkg}_INCLUDE_DIR)
    list(GET ${pkg}_INCLUDE_DIR 0 TMP_DIR)
    get_filename_component(TMP_DIR ${TMP_DIR} PATH)
    get_filename_component(TMP_DIR ${TMP_DIR} PATH)
  endif()
  if(NOT TMP_DIR OR NOT EXISTS ${TMP_DIR})
    message("${pkg}_INCLUDE_DIR is not set")
  else()
    set(SV_${pkg}_DIR ${TMP_DIR})
  endif()
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# System Macros
macro(env_variable_to_value_variable value_variable variable)
	if(WIN32 AND NOT UNIX)
		set(${value_variable} "%${variable}%")
	endif()
	if(UNIX)
		set(${value_variable} "$${variable}")
	endif()
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
function(append_env_string evn_var value output_variable)
	env_variable_to_value_variable(ENV_VALUE ${evn_var})
	set(${output_variable} "${ENV_SET_COMMAND} ${evn_var}=${ENV_VALUE}${ENV_SEPERATOR}${value}" PARENT_SCOPE)
endfunction()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
function(set_env_string evn_var value output_variable)
	set(${output_variable} "${ENV_SET_COMMAND} ${evn_var}=${value}\n" PARENT_SCOPE)
endfunction()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
macro(set_env_string_concat evn_var value output_variable)
	set_env_string(${evn_var} ${value} _tmp)
	set(${output_variable} "${${output_variable}}${_tmp}")
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
macro(append_env_string_concat evn_var value output_variable)
	append_env_string(${evn_var} ${value} _tmp)
	set(${output_variable} "${${output_variable}}${_tmp}\n")
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
macro(simvascular_property_list_find_and_replace TARGET PROPERTY VALUE NEWVALUE)
  get_target_property(_LIST_VAR ${TARGET} ${PROPERTY})
  if("${_LIST_VAR}" STREQUAL "_LIST_VAR-NOTFOUND")
    dev_message("Property list find and replace, property ${PROPERTY} not found")
  else()
    simvascular_list_find_and_replace(_LIST_VAR "${VALUE}" ${NEWVALUE})
    set_target_properties(${TARGET} PROPERTIES ${PROPERTY} "${_LIST_VAR}")
  endif()
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
macro(simvascular_list_find_and_replace LIST VALUE NEWVALUE)
  set(_COUNT "0")
  foreach(ITEM ${${LIST}})
    string(REGEX MATCH "${VALUE}" _FOUND ${ITEM})
    if(_FOUND)
      simvascular_list_replace(${LIST} ${_COUNT} ${NEWVALUE})
    endif()
    math(EXPR _COUNT "${_COUNT}+1")
  endforeach()
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
macro(simvascular_list_replace LIST INDEX NEWVALUE)
  list(INSERT ${LIST} ${INDEX} ${NEWVALUE})
  math(EXPR __INDEX "${INDEX} + 1")
  list (REMOVE_AT ${LIST} ${__INDEX})
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
function(print_target_properties tgt)
  if(NOT TARGET ${tgt})
    message("There is no target named '${tgt}'")
      return()
    endif()

    execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)
    # Convert command output into a CMake list
    string(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
    string(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")

    foreach (prop ${CMAKE_PROPERTY_LIST})
      string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" prop ${prop})
      # message ("Checking ${prop}")
      get_property(propval TARGET ${tgt} PROPERTY ${prop} SET)
      if (propval)
        get_target_property(propval ${tgt} ${prop})
        message ("${tgt} ${prop} = ${propval}")
      endif()
   endforeach(prop)
endfunction(print_target_properties)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
macro(simvascular_get_major_minor_version version major_version minor_version)
  string(REPLACE "." ";" version_list ${version})
  list(GET version_list 0 ${major_version})
  list(GET version_list 1 ${minor_version})
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
macro(simvascular_today YEAR MONTH DAY)
  if(WIN32)
    execute_process(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE RESULT)
    string(REGEX REPLACE "(..)/(..)/(....).*" "\\3;\\2;\\1" RESULT ${RESULT})
  elseif(UNIX)
    execute_process(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE RESULT)
    string(REGEX REPLACE "(..)/(..)/(....).*" "\\3;\\2;\\1" RESULT ${RESULT})
  else(WIN32)
    message(SEND_ERROR "date not implemented")
    set(YEAR 0000)
    set(MONTH 00)
    set(DAY 00)
  endif (WIN32)
  list(GET RESULT 0 YEAR)
  list(GET RESULT 1 MONTH)
  list(GET RESULT 2 DAY)
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# simvascular_add_new_external
macro(simvascular_add_new_external proj version use shared dirname)
  option(SV_USE_${proj} "Enable ${proj} Plugin" ${use})
  option(SV_USE_${proj}_SHARED "Build ${proj} libraries as shared libs" ${shared})
  option(SV_USE_${proj}_QT_GUI "Build ${proj} in the Qt GUI" ${use})
  mark_as_advanced(SV_USE_${proj}_QT_GUI)

  set(${proj}_VERSION "${version}" CACHE TYPE STRING)
  simvascular_get_major_minor_version(${${proj}_VERSION} ${proj}_MAJOR_VERSION ${proj}_MINOR_VERSION)
  set(SV_EXT_${proj}_SRC_DIR ${SV_EXTERNALS_SRC_DIR}/${dirname}-${${proj}_VERSION})
  set(SV_EXT_${proj}_BIN_DIR ${SV_EXTERNALS_BIN_DIR}/${dirname}-${${proj}_VERSION})
  set(SV_EXT_${proj}_BLD_DIR ${SV_EXTERNALS_BLD_DIR}/${dirname}-${${proj}_VERSION})
  set(SV_EXT_${proj}_PFX_DIR ${SV_EXTERNALS_PFX_DIR}/${dirname}-${${proj}_VERSION})

  # Install rules
  set(SV_EXTERNALS_${proj}_INSTALL_PREFIX ${SV_EXTERNALS_INSTALL_PREFIX}/${dirname}-${${proj}_VERSION})
  if(${CMAKE_PROJECT_NAME}_ENABLE_DISTRIBUTION)
    set(LIB_DESTINATION "${SV_EXTERNALS_INSTALL_PREFIX}")
  else()
    set(LIB_DESTINATION "${SV_EXTERNALS_${proj}_INSTALL_PREFIX}")
  endif()

  if(NOT SV_INSTALL_${proj}_RUNTIME_DIR)
    set(SV_INSTALL_${proj}_RUNTIME_DIR ${LIB_DESTINATION}/bin)
  endif()

  if(NOT SV_INSTALL_${proj}_LIBRARY_DIR)
    set(SV_INSTALL_${proj}_LIBRARY_DIR ${LIB_DESTINATION}/lib)
  endif()

  if(NOT SV_INSTALL_${proj}_ARCHIVE_DIR)
    set(SV_INSTALL_${proj}_ARCHIVE_DIR ${LIB_DESTINATION}/lib)
  endif()

  if(NOT SV_INSTALL_${proj}_INCLUDE_DIR)
    set(SV_INSTALL_${proj}_INCLUDE_DIR ${LIB_DESTINATION}/include)
  endif()

  if(SV_USE_${proj})
    list(APPEND SV_EXTERNALS_LIST ${proj})
    set(SV_${proj}_DIR ${SV_EXT_${proj}_BIN_DIR})
    if(NOT ${proj}_DIR)
      set(${proj}_DIR "" CACHE PATH "For external projects with a Config.cmake file, path to that file; for externals without a Config.cmake, the path to the toplevel bin directory")
    else()
      set(${proj}_DIR "${${proj}_DIR}" CACHE PATH "For external projects with a Config.cmake file, path to that file; for externals without a Config.cmake, the path to the toplevel bin directory")
    endif()
    if("${proj}" STREQUAL "MITK")
      if(SV_USE_MITK_CONFIG)
        set(SV_${proj}_DIR ${SV_EXT_${proj}_BLD_DIR}/MITK-build)
      endif()
    endif()
  endif()
  if(SV_EXTERNALS_LIST)
    list(REMOVE_DUPLICATES SV_EXTERNALS_LIST)
  endif()
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# sv_externals_add_new_external
# \brief Create new external and set variables with default values based on inputs
macro(sv_externals_add_new_external proj version use shared dirname install_dirname)
  option(SV_EXTERNALS_ENABLE_${proj} "Enable ${proj} Plugin" ${use})
  option(SV_EXTERNALS_ENABLE_${proj}_SHARED "Build ${proj} libraries as shared libs" ${shared})
  mark_as_advanced(SV_EXTERNALS_ENABLE_${proj}_SHARED)
  option(SV_EXTERNALS_DOWNLOAD_${proj} "Download instead of build ${proj}; Unused for tcl, tk, tcllib, tklib, numpy, pip" ${use})
  mark_as_advanced(SV_EXTERNALS_DOWNLOAD_${proj})

  # Version
  set(SV_EXTERNALS_${proj}_VERSION "${version}" CACHE TYPE STRING)
  mark_as_advanced(SV_EXTERNALS_${proj}_VERSION)
  simvascular_get_major_minor_version(${SV_EXTERNALS_${proj}_VERSION} SV_EXTERNALS_${proj}_MAJOR_VERSION SV_EXTERNALS_${proj}_MINOR_VERSION)

  # Src, bin, build, prefic dirs
  set(${proj}_VERSION_DIR ${dirname}-${SV_EXTERNALS_${proj}_VERSION})
  set(SV_EXTERNALS_${proj}_SRC_DIR ${SV_EXTERNALS_TOPLEVEL_SRC_DIR}/${${proj}_VERSION_DIR})
  if(NOT "${install_dirname}" STREQUAL "none")
    set(SV_EXTERNALS_${proj}_BIN_DIR ${SV_EXTERNALS_TOPLEVEL_BIN_DIR}/${install_dirname}-${SV_EXTERNALS_${proj}_VERSION})
  else()
    set(SV_EXTERNALS_${proj}_BIN_DIR ${SV_EXTERNALS_TOPLEVEL_BIN_DIR}/${${proj}_VERSION_DIR})
  endif()
  set(SV_EXTERNALS_${proj}_BLD_DIR ${SV_EXTERNALS_TOPLEVEL_BLD_DIR}/${${proj}_VERSION_DIR})
  set(SV_EXTERNALS_${proj}_PFX_DIR ${SV_EXTERNALS_TOPLEVEL_PFX_DIR}/${${proj}_VERSION_DIR})

  # Install dirs
  if(NOT SV_EXTERNALS_INSTALL_${proj}_RUNTIME_DIR)
    set(SV_EXTERNALS_INSTALL_${proj}_RUNTIME_DIR ${SV_EXTERNALS_${proj}_BIN_DIR}/bin)
  endif()

  if(NOT SV_EXTERNALS_INSTALL_${proj}_LIBRARY_DIR)
    set(SV_EXTERNALS_INSTALL_${proj}_LIBRARY_DIR ${SV_EXTERNALS_${proj}_BIN_DIR}/lib)
  endif()

  if(NOT SV_EXTERNALS_INSTALL_${proj}_ARCHIVE_DIR)
    set(SV_EXTERNALS_INSTALL_${proj}_ARCHIVE_DIR ${SV_EXTERNALS_${proj}_BIN_DIR}/lib)
  endif()

  if(NOT SV_EXTERNALS_INSTALL_${proj}_INCLUDE_DIR)
    set(SV_EXTERNALS_INSTALL_${proj}_INCLUDE_DIR ${SV_EXTERNALS_${proj}_BIN_DIR}/include)
  endif()

  # Add to externals list
  if(SV_EXTERNALS_ENABLE_${proj})
    list(APPEND SV_EXTERNALS_LIST ${proj})
  endif()
  if(SV_EXTERNALS_LIST)
    list(REMOVE_DUPLICATES SV_EXTERNALS_LIST)
  endif()

  # Add install step for each external
  if(NOT "${install_dirname}" STREQUAL "none")
    simvascular_today(YEAR MONTH DAY)
    set(SV_EXTERNALS_${proj}_TAR_INSTALL_NAME ${SV_PLATFORM_DIR}.${SV_PLATFORM_VERSION_DIR}.${SV_COMPILER_DIR}.${SV_COMPILER_VERSION_DIR}.${SV_ARCH_DIR}.${SV_BUILD_TYPE_DIR}.${YEAR}.${MONTH}.${DAY}.${install_dirname}.${SV_EXTERNALS_${proj}_VERSION})
    if(EXISTS "${SV_EXTERNALS_TAR_INSTALL_DIR}")
      install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E tar -czvf ${SV_EXTERNALS_TAR_INSTALL_DIR}/${SV_EXTERNALS_${proj}_TAR_INSTALL_NAME}.tar.gz ${SV_EXTERNALS_${proj}_BIN_DIR}
        WORKING_DIRECTORY ${SV_EXTERNALS_TOPLEVEL_BIN_DIR})")
    endif()
  endif()

  # Set up download stuff if downloading
  if(NOT "${install_dirname}" STREQUAL "none")
    if(SV_EXTERNALS_DOWNLOAD_${proj})
      set(${proj}_TEST_FILE "${SV_EXTERNALS_URL}/${SV_KERNEL_DIR}/${SV_PLATFORM_DIR}/externals_compiler_info.txt")
      file(DOWNLOAD "${${proj}_TEST_FILE}" "${SV_EXTERNALS_${proj}_PFX_DIR}/externals_compiler_info.txt" STATUS _status LOG _log INACTIVITY_TIMEOUT 5 TIMEOUT 5)
      list(GET _status 0 err)
      list(GET _status 1 msg)
      if(err)
        message(FATAL_ERROR "The operating system does not have any available per-built binaries. See the build documentation to build your own.")
      else()
        simvascular_read_file("${SV_EXTERNALS_${proj}_PFX_DIR}/externals_compiler_info.txt" FILE_CONTENTS)
        sv_externals_check_versioning("${FILE_CONTENTS}" ${SV_PLATFORM_VERSION_DIR} ${SV_COMPILER_DIR} ${SV_COMPILER_VERSION_DIR} SV_DOWNLOAD_DIR)
        string(REPLACE "/" "." SV_TAR_PREFIX "${SV_DOWNLOAD_DIR}")
        set(SV_EXTERNALS_${proj}_BINARIES_URL "${SV_EXTERNALS_URL}/${SV_KERNEL_DIR}/${SV_PLATFORM_DIR}/${SV_DOWNLOAD_DIR}/${SV_PLATFORM_DIR}.${SV_TAR_PREFIX}.${install_dirname}.${SV_EXTERNALS_${proj}_VERSION}.tar.gz")
      endif()
    endif()
  endif()
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# sv_externals_read_file
macro(simvascular_read_file file_name file_contents)
  file(READ "${file_name}" ${file_contents})
  # Convert file contents into a CMake list (where each element in the list
  # is one line of the file)
  #
  string(REGEX REPLACE ";" "\\\\;" ${file_contents} "${${file_contents}}")
  string(REGEX REPLACE "\n" ";" ${file_contents} "${${file_contents}}")
endmacro()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# sv_externals_check_versioning
macro(sv_externals_check_versioning check_file_contents platform_version compiler compiler_version output_dir)

  # Initiate loop variables
  set(found_options complete_match
                    platform_compiler
                    platform_only
                    compiler_only
                    compiler_type_only
                    nothing)
  foreach(foption ${found_options})
    set(${foption} FALSE)
    set(${foption}_oldest_platform_ver "100000000")
    set(${foption}_oldest_compiler_ver "100000000")
    set(${foption}_oldest_compiler "")
    set(${foption}_rest_of_line "")
  endforeach()

  # Loop through file contents
  foreach(line ${check_file_contents})

    # Seperate out the platform version, compiler dir, and the date
    string(REPLACE "/" ";" line_list ${line})
    list(GET line_list 0 fileline_platform_version)
    list(GET line_list 1 fileline_compiler)
    list(GET line_list 2 fileline_compiler_version)
    list(GET line_list 3 4 5 fileline_rest_of_line_list)
    string(REPLACE ";" "/" fileline_rest_of_line "${fileline_rest_of_line_list}")

    # Check what matches
    if("${platform_version}/${compiler}/${compiler_version}" STREQUAL "${fileline_platform_version}/${fileline_compiler}/${fileline_compiler_version}")

      # Exact match!
      set(complete_match TRUE)
      set(complete_match_rest_of_line "${fileline_rest_of_line}")
      break()

    elseif("${platform_version}/${compiler}" STREQUAL "${fileline_platform_version}/${fileline_compiler}")

      # Platform and compiler match!
      set(platform_compiler TRUE)
      if(fileline_compiler_version VERSION_LESS "${platform_compiler_oldest_compiler_ver}")
        set(platform_compiler_oldest_compiler_ver "${fileline_compiler_version}")
        set(platform_compiler_rest_of_line "${fileline_rest_of_line}")
      endif()

    elseif("${compiler}/${compiler_version}" STREQUAL "${fileline_compiler}/${fileline_compiler_version}")

      # Just the compiler dir matches, this is actually pretty good
      set(compiler_only TRUE)
      if(fileline_platform_version VERSION_LESS "${compiler_only_oldest_platform_ver}")
        set(compiler_only_oldest_platform_ver "${fileline_platform_version}")
        set(compiler_only_rest_of_line "${fileline_rest_of_line}")
      endif()

    elseif("${compiler}" STREQUAL "${fileline_compiler}")

      # Just the compiler type was found, ehh okay
      set(compiler_type_only TRUE)
      if(fileline_compiler_version VERSION_LESS "${compiler_type_only_oldest_compiler_ver}")
        set(compiler_type_only_oldest_compiler_ver "${fileline_compiler_version}")
        set(compiler_type_only_oldest_platform_ver "${fileline_platform_version}")
        set(compiler_type_only_rest_of_line "${fileline_rest_of_line}")
      endif()

    elseif("${platform_version}" STREQUAL "${fileline_platform_version}")

      # Only platform found, use an arbitrary compiler
      set(platform_only TRUE)
      if(fileline_compiler_version VERSION_LESS "${platform_only_oldest_compiler_ver}")
        set(platform_only_oldest_compiler_ver "${fileline_compiler_version}")
        set(platform_only_oldest_compiler "${fileline_compiler}")
        set(platform_only_rest_of_line "${fileline_rest_of_line}")
      endif()

    else()

      # Nothing found this line, set oldest platform version
      set(nothing TRUE)
      if("${fileline_platform_version}" VERSION_LESS "${nothing_oldest_platform_ver}")
        set(nothing_oldest_platform_ver "${fileline_platform_version}")
        set(nothing_oldest_compiler "${fileline_compiler}")
        set(nothing_oldest_compiler_ver "${fileline_compiler_version}")
        set(nothing_rest_of_line "${fileline_rest_of_line}")
      endif()

    endif()

  endforeach()

  # Set the generic warning
  set(GENERIC_MESSAGE "Pre-built binaries for the operating system and compiler do not exist! The best possible match will be downloaded; however, problems may occur, especially if the pre-built binaries are compiled with a different compiler.")

  # Find what happened in loop!
  if(complete_match)
    # Simple, everything matches, yeah!
    set(${output_dir} "${platform_version}/${compiler}/${compiler_version}/${complete_match_rest_of_line}")

  elseif(platform_compiler)
    # Compiler found but wrong version
    message(WARNING "${GENERIC_MESSAGE} Pre-built binaries for ${SV_PLATFORM_DIR} version ${platform_version} and compiler ${compiler}/${platform_compiler_oldest_compiler_ver} are being downloaded and used. Proceed with caution!")
    set(${output_dir} "${platform_version}/${compiler}/${platform_compiler_oldest_compiler_ver}/${platform_compiler_rest_of_line}")

  elseif(compiler_only)
    # Platform not there, but compiler was!
    message(WARNING "${GENERIC_MESSAGE} Pre-built binaries for ${SV_PLATFORM_DIR} version ${compiler_only_oldest_platform_ver} and compiler ${compiler}/${compiler_version} are being downloaded and used. Proceed with caution!")
    set(${output_dir} "${compiler_only_oldest_platform_ver}/${compiler}/${compiler_version}/${compiler_only_rest_of_line}")

  elseif(compiler_type_only)
    # Platform not there, but compiler was!
    message(WARNING "${GENERIC_MESSAGE} Pre-built binaries for ${SV_PLATFORM_DIR} version ${compiler_type_only_oldest_platform_ver} and compiler ${compiler}-${compiler_type_only_oldest_compiler_ver} are being downloaded and used. Proceed with caution!")
    set(${output_dir} "${compiler_type_only_oldest_platform_ver}/${compiler}/${compiler_type_only_oldest_compiler_ver}/${compiler_type_only_rest_of_line}")

  elseif(platform_only)
    # Even worse, issue warning and leave
    message(WARNING "${GENERIC_MESSAGE} Pre-built binaries for ${SV_PLATFORM_DIR} version ${platform_version} and compiler ${platform_only_oldest_compiler}-${platform_only_oldest_compiler_ver} are being downloaded and used. Proceed with caution!")
    set(${output_dir} "${platform_version}/${platform_only_oldest_compiler}/${platform_only_oldest_compiler_ver}/${platform_only_rest_of_line}")
  else()
    # The worst! fatal error
    message(WARNING "${GENERIC_MESSAGE} Pre-built binaries for ${SV_PLATFORM_DIR} version ${nothing_oldest_platform_ver} and compiler ${nothing_oldest_compiler}/${nothing_oldest_compiler_ver}")
    set(${output_dir} "${nothing_oldest_platform_ver}/${nothing_oldest_compiler}-${nothing_oldest_compiler_ver}/${nothing_rest_of_line}")
  endif()

endmacro()
#-----------------------------------------------------------------------------
