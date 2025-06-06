# Create the Lorenz 63 Fortran executable
add_executable(lorenz63 main.f90)

# Link to required libraries
target_link_libraries(lorenz63
    PRIVATE
        metada::backends::common::utils
        lorenz63_state
        lorenz63_geometry
        lorenz63_model)

# Make include directories available
target_include_directories(lorenz63
    PRIVATE
        ${CMAKE_BINARY_DIR}/modules)

# Create the C++ wrapper example executable
add_executable(lorenz63_cpp_example lorenz63_cpp_example.cpp)

# Link to required libraries
target_link_libraries(lorenz63_cpp_example
    PRIVATE
        metada::traits
        metada::framework::runs
        metada::framework::interfaces
        metada::framework::adapters
        metada::backends::common::utils
        lorenz63_state
        lorenz63_geometry
        lorenz63_model)

# Make include directories available
target_include_directories(lorenz63_cpp_example
    PRIVATE
        ${CMAKE_BINARY_DIR}/modules
        ${CMAKE_SOURCE_DIR}/src/backends/interfaces)

AddFormatTarget(lorenz63_cpp_example ${CMAKE_CURRENT_SOURCE_DIR})