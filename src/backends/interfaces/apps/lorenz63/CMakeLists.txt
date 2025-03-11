# Create the Lorenz 63 Fortran executable
add_executable(lorenz63 main.f90)

# Link to required libraries
target_link_libraries(lorenz63 
  PRIVATE 
    lorenz63_state
    lorenz63_geometry
    lorenz63_model)

# Make include directories available
target_include_directories(lorenz63 
  PRIVATE
    ${CMAKE_BINARY_DIR}/modules)

# Add C++ example if C++ is enabled
if(CMAKE_CXX_COMPILER)
  # Create the C++ wrapper example executable
  add_executable(lorenz63_cpp_example lorenz63_cpp_example.cpp)
  
  # Link to required libraries
  target_link_libraries(lorenz63_cpp_example
    PRIVATE
      metada::framework::runs
      metada::framework::traits
      metada::framework::common
      metada::backends::common::utils::logger
      metada::backends::common::utils::config
      lorenz63_state
      lorenz63_geometry
      lorenz63_model)
  
  # Make include directories available
  target_include_directories(lorenz63_cpp_example
    PRIVATE
      ${CMAKE_BINARY_DIR}/modules
      ${CMAKE_SOURCE_DIR}/src/backends/interfaces)
  
endif()