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


#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.3.0)

if (POLICY CMP0040)
  cmake_policy(SET CMP0040 OLD)
endif()
if (POLICY CMP0046)
  cmake_policy(SET CMP0046 OLD)
endif()
if (POLICY CMP0053)
  cmake_policy(SET CMP0053 OLD)
endif()
if (POLICY CMP0054)
  cmake_policy(SET CMP0054 OLD)
endif()
if (POLICY CMP0025) # Compiler ID when using clang set to clang not appleclang
  cmake_policy(SET CMP0025 OLD)
endif()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
Project(SV)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake"
  "${CMAKE_MODULE_PATH}")
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# CMake Includes
include(CheckLibraryExists)
include(GetPrerequisites)
include(GenerateExportHeader)
include(ExternalProject)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# SimVascular Includes
include(SimVascularMacros)
include(SimVascularFunctions)
include(SimVascularSystemSetup)
include(SimVascularFunctionCheckCompilerFlags)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# RPATH handling
# Skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  TRUE)
# Use the install RPATH
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# The RPATH that will be used when installing
set(CMAKE_INSTALL_RPATH "")
# Automatically add determined parts of the RPATH which point to
# directories outside the build tree to the install tree
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Include Initial Options
include(SimVascularOptions)
include(SimVascularDependentOptions)
include(SimVascularInstallSetup)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Set EXE names for the executables
set(SV_PRESOLVER_EXE "svpre")
set(SV_FLOWSOLVER_EXE "svsolver")
set(SV_POSTSOLVER_EXE "svpost")
set(SV_SVFSI_EXE "svFSI")
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Include SimVascular Internals
include(SimVascularInternals)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Externals
include(SimVascularExternals)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# ThirdParty
include(SimVascularThirdParty)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Add GLOBAL_DEFINES to definitions
dev_message(" Global Flags: ${GLOBAL_DEFINES}")
if(GLOBAL_DEFINES)
  list(REMOVE_DUPLICATES GLOBAL_DEFINES)
endif()
add_definitions(${GLOBAL_DEFINES})
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Add Local Includes
set(LOCAL_LIBDIR Lib)
set(SV_SCRIPT_TARGETS "" CACHE INTERNAL "" FORCE)
set(INDENTS "")
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Include the Core subdirs
# These are handled here because other sub directories depend on them
# set
set(SV_INCLUDE_DIRS "")
set(SV_INCLUDE_DIRS ${SV_INCLUDE_DIRS}
  "FlowSolvers/Include"
  )

foreach(dir ${SV_INCLUDE_DIRS})
  add_subdirectory(${dir})
  include_directories(${SV_SOURCE_DIR}/${dir} ${SV_BINARY_DIR}/${dir})
endforeach()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Add directories
# Modules
add_subdirectory(FlowSolvers)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Enable distribution if turned on
if(SV_ENABLE_DISTRIBUTION)
  add_subdirectory(${SV_DISTRIBUTION_DIR} ${SV_BINARY_DISTRIBUTION_DIR})
endif()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Configure Exe Scripts, this should be the last subdirectory
include(SimVascularInstallLibs)
add_subdirectory(Scripts)
#-----------------------------------------------------------------------------
