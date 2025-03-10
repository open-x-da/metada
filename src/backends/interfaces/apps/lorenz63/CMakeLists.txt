# Create the Lorenz 63 executable
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

# Install the executable to the correct location
install(TARGETS lorenz63
  RUNTIME DESTINATION ${METADA_INSTALL_BINDIR})