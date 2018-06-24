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

macro(simvascular_install_prereqs tar location component)
  get_prerequisites(${tar} DEPENDENCIES 1 1 "" "")
  foreach(DEPENDENCY ${DEPENDENCIES})
    get_filename_component(DEPENDENCY_NAME "${DEPENDENCY}" NAME)
    get_filename_component(DEPENDENCY_NAME_WE "${DEPENDENCY}" NAME_WE)
    get_filename_component(DEPENDENCY_ACTUAL "${DEPENDENCY}" REALPATH)
    dev_message("[${tar} Dep] ${DEPENDENCY_ACTUAL}")
    install(FILES "${DEPENDENCY_ACTUAL}"
      DESTINATION "${location}"
      RENAME "${DEPENDENCY_NAME}"
      COMPONENT "${component}"
      )
  endforeach()
endmacro()

foreach(var ${SV_EXTERNAL_SHARED_LIBS})
  if(NOT ${var}_DLL_LIBRARIES)
    set(${var}_DLL_LIBRARIES ${${var}_LIBRARY} ${${var}_LIBRARIES} ${${var}_DLL_LIBRARY})
  endif()
  if(${var}_DLL_LIBRARIES)
    list(REMOVE_DUPLICATES ${var}_DLL_LIBRARIES)
  endif()
  if(SV_ENABLE_DISTRIBUTION OR NOT SV_USE_SYSTEM_${var} )
    # We only want to install libraries they are not system libraries.
    # we want to install them regardles if this is a distribution.
    dev_message("[${var}]  This is either not a system library, or this is a distribution install.")
    if(NOT SV_INSTALL_${var}_RUNTIME_DIR)
      dev_message("[${var}]  SV_INSTALL_${var}_RUNTIME_DIR does not exist, ${lib_name} installing into externals")
      # Variables with a install directory are handled differently
      # if the directory doesn't exist, install in in the catch all
      foreach(lib ${${var}_DLL_LIBRARIES})
        dev_message("[${var}]  Do we install ${lib}")
        get_filename_component(lib_name ${lib} NAME)
        get_filename_component(_EXT ${lib_name} EXT)
        get_filename_component(lib_name_we ${lib} NAME_WE)
        get_filename_component(lib_path ${lib} PATH)
        if(_EXT MATCHES "${CMAKE_SHARED_LIBRARY_SUFFIX}$")
          dev_message("[${var}] ${lib_name} marked for install into ${SV_INSTALL_EXTERNALS_RUNTIME_DIR}")
          if(IS_SYMLINK ${lib})
            dev_message("${lib_name} ${lib_name_we} is a symbolic link\n${lib_path}/${lib_name_we}*")
            #get_filename_component(lib ${})
            FILE(GLOB lib "${lib_path}/${lib_name_we}*")
          endif()
          dev_message("Installing ${lib}")
          install(FILES ${lib} DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT ExternalLibraries)
          simvascular_install_prereqs(${lib} ${SV_INSTALL_EXTERNALS_RUNTIME_DIR} ExternalLibraries)
        else()
          dev_message("[${var}]  ${lib_name} is not a shared lib, removing ${lib_name}, EXT ${_EXT}, do not install")
          list(REMOVE_ITEM ${var}_DLL_LIBRARIES ${lib})
        endif()
      endforeach()
    else()
      #debug only
      dev_message("SV_INSTALL_${var}_RUNTIME_DIR exists, ${var} will be handled elsewhere")
    endif()
  else()
    dev_message("[${var}] This is system library, and distribution is not enabled.")
  endif()
  dev_message("")
endforeach()

#-----------------------------------------------------------------------------
# MPI
if(SV_USE_MPI AND NOT SV_USE_DUMMY_MPI)
  if(SV_ENABLE_DISTRIBUTION OR NOT SV_USE_SYSTEM_MPI)
    if(NOT WIN32)
      # MPI needs to be insalled, if its not a system LIBRARY.
      # or this is a distribution, unless its WIN32

      # Find all the libraries and install them
      getListOfVars("MPI" "LIBRARIES" _VARLIST)
      getListOfVars_concat("MPI" "LIBRARY" _VARLIST)
      listvars2vals(_VARLIST MPI_LIBS)
      set(MPI_INSTALL_LIBS)
      foreach(lib ${MPI_LIBS})
        get_filename_component(lib_name ${lib} NAME)
        get_filename_component(_EXT ${lib_name} EXT)
        get_filename_component(lib_name_we ${lib} NAME_WE)
        get_filename_component(lib_path ${lib} PATH)
        get_filename_component(real_lib ${lib} REALPATH)
        file(GLOB libs "${lib_path}/${lib_name_we}*")
        foreach(lib2 ${libs})
          if(NOT ${lib2} MATCHES "${CMAKE_STATIC_LIBRARY_SUFFIX}$")
            set(MPI_INSTALL_LIBS ${MPI_INSTALL_LIBS} ${lib2})
          endif()
        endforeach()
      endforeach()
      if(MPI_INSTALL_LIBS)
        list(REMOVE_DUPLICATES MPI_INSTALL_LIBS)
      endif()
      install(FILES ${MPI_INSTALL_LIBS}
        DESTINATION ${SV_INSTALL_MPI_LIBRARY_DIR} COMPONENT MPILibraries)
      foreach(lib ${MPI_INSTALL_LIBS})
        simvascular_install_prereqs(${lib} ${SV_INSTALL_MPI_RUNTIME_DIR} MPILibraries)
    endforeach()
    endif()
    #find MPIEXEC's path, and install it and everything with MPI or hydra in the name
    get_filename_component(MPIEXEC_PATH ${MPIEXEC} PATH)
    file(GLOB MPIEXEC_INSTALL_PROGS "${MPIEXEC_PATH}/mpiexec*")
    file(GLOB HYDRA_INSTALL_PROGS "${MPIEXEC_PATH}/*hydra*")
    set(MPIEXEC_INSTALL_PROGS ${MPIEXEC_INSTALL_PROGS} ${HYDRA_INSTALL_PROGS} ${MPIEXEC})
    if(MPIEXEC_INSTALL_PROGS)
      list(REMOVE_DUPLICATES MPIEXEC_INSTALL_PROGS)
      # CHANGE NO TARGETS
      # install(PROGRAMS ${MPIEXEC_INSTALL_PROGS} DESTINATION ${SV_INSTALL_MPI_EXE_DIR})
      install(PROGRAMS ${MPIEXEC_INSTALL_PROGS}
        DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT MPIExecutables)
    endif()
    if(APPLE)
      install(CODE "
      file(GLOB INSTALLED_MPIEXECS \"\${CMAKE_INSTALL_PREFIX}/${SV_INSTALL_RUNTIME_DIR}/mpiexec-*\")
      if(INSTALLED_MPIEXECS)
        list(GET INSTALLED_MPIEXECS 0 WORK_MPIEXEC)
        get_filename_component(ONE_MPIEXEC \"\${WORK_MPIEXEC}\" NAME)
        if (EXISTS \"\${CMAKE_INSTALL_PREFIX}/${SV_INSTALL_RUNTIME_DIR}/mpiexec\")
          file(REMOVE \"\${CMAKE_INSTALL_PREFIX}/${SV_INSTALL_RUNTIME_DIR}/mpiexec\")
        endif()
        execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink \"\${ONE_MPIEXEC}\" \"mpiexec\" WORKING_DIRECTORY \${CMAKE_INSTALL_PREFIX}/${SV_INSTALL_RUNTIME_DIR})
      endif()
      "
      COMPONENT MPIExecutables)
    endif()
  endif()
endif()

