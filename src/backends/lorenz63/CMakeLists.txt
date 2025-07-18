cmake_minimum_required(VERSION 3.10)

# Add subdirectories for Lorenz 63 components
add_subdirectory(state)
add_subdirectory(geometry)
add_subdirectory(model)
add_subdirectory(apps)

# Create interfaces backend library
add_library(interfaces_backend INTERFACE)

# Add alias for library
add_library(interfaces::backend
    ALIAS interfaces_backend)

# Set include directories
target_include_directories(interfaces_backend
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})

# Link dependencies
target_link_libraries(interfaces_backend
    INTERFACE
        metada::framework::interfaces
        metada::framework::tools::config
        lorenz63_state
        lorenz63_geometry
        lorenz63_model)

# Add target to build all Lorenz 63 components (for convenience)
add_custom_target(lorenz63_all
  DEPENDS
  lorenz63_state
  lorenz63_geometry
  lorenz63_model
  lorenz63)