# Set Coronary Flag
if(SV_USE_CORONARY)
  set(VER_CORONARY 1)
else()
  set(VER_CORONARY 0)
endif()
# Set Closed Loop Flag
if(SV_USE_CLOSEDLOOP)
  set(VER_CLOSEDLOOP 1)
else()
  set(VER_CLOSEDLOOP 0)
endif()
# Set Variable wall Flag
if(SV_USE_VARWALL)
  set(VER_VARWALL 1)
else()
  set(VER_VARWALL 0)
endif()
# Set Use VTK Flag
if(SV_USE_VTK)
  set(VER_USE_VTK 1)
else()
  set(VER_USE_VTK 0)
endif()

configure_file(simvascular_version.h.in ${CMAKE_CURRENT_BINARY_DIR}/simvascular_version.h)
configure_file(simvascular_options.h.in ${CMAKE_CURRENT_BINARY_DIR}/simvascular_options.h)

# Write include file
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cvFlowsolverOptions.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cvFlowsolverOptions.h"
)
