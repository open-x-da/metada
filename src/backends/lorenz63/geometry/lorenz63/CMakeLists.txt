# Create the Lorenz 63 geometry library
add_library(lorenz63_geometry STATIC 
  geometry.f90
  geometry_c_bindings.f90)

# Set properties for the target
set_target_properties(lorenz63_geometry PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Link to required libraries
target_link_libraries(lorenz63_geometry 
  PUBLIC 
    lorenz63_state)

# Make include directories available
target_include_directories(lorenz63_geometry 
  PUBLIC
    ${CMAKE_BINARY_DIR}/modules)

# Add library alias
add_library(metada::backends::lorenz63::geometry
  ALIAS lorenz63_geometry)

AddFormatTarget(lorenz63_geometry ${CMAKE_CURRENT_SOURCE_DIR})