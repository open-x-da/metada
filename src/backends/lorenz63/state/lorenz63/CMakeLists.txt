# Create the Lorenz 63 state library
add_library(lorenz63_state STATIC 
  state.f90
  state_c_bindings.f90)

# Set properties for the target
set_target_properties(lorenz63_state PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Make include directories available to dependent targets
target_include_directories(lorenz63_state 
  PUBLIC
    ${CMAKE_BINARY_DIR}/modules)

# Add library alias
add_library(metada::backends::lorenz63::state
  ALIAS lorenz63_state)

AddFormatTarget(lorenz63_state ${CMAKE_CURRENT_SOURCE_DIR})
AddFortranFormatTarget(lorenz63_state ${CMAKE_CURRENT_SOURCE_DIR})