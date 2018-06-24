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
# getListOfVars -
#
function(getListOfVars _prefix _suffix _varResult)
  get_cmake_property(_vars VARIABLES)
  string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*${_suffix}" _matchedVars "${_vars}")
	set (${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

function(getListOfVars_concat _prefix _suffix _varResult)
  get_cmake_property(_vars VARIABLES)
  string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*${_suffix}" _matchedVars "${_vars}")
  set (${_varResult} ${${_varResult}} ${_matchedVars} PARENT_SCOPE)
endfunction()

function(getListOfVarsPrefix _prefix _varResult)
  get_cmake_property(_vars VARIABLES)
  string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
  set (${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

function(getListofVarsCat _varResult)
  set(options ) 
  set(oneValueArgs)
  set(multiValueArgs SUFFIXES PREFIXES)
  CMAKE_PARSE_ARGUMENTS("" 
    "${options}"
    "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  set(_ACCUM_VARLIST)
  foreach(pre ${_PREFIXES})
    foreach(suf ${_SUFFIXES})
      getListOfVars("${pre}" "${suf}" _VARLIST)
      set(_ACCUM_VARLIST ${_ACCUM_VARLIST} ${_VARLIST})
    endforeach()
  endforeach()
  set(_varResult ${_ACCUM_VARLIST} PARENT_SCOPE)
endfunction()

