# Create the Lorenz 63 model library
add_library(lorenz63_model STATIC 
  model.f90
  model_c_bindings.f90)

# Set properties for the target
set_target_properties(lorenz63_model PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Link to required libraries
target_link_libraries(lorenz63_model 
  PUBLIC 
    lorenz63_state)

# Make include directories available
target_include_directories(lorenz63_model 
  PUBLIC
    ${CMAKE_BINARY_DIR}/modules)

# Add library alias
add_library(metada::backends::lorenz63::model
  ALIAS lorenz63_model)

AddFormatTarget(lorenz63_model ${CMAKE_CURRENT_SOURCE_DIR})
AddFortranFormatTarget(lorenz63_model ${CMAKE_CURRENT_SOURCE_DIR})