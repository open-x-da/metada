# Add LETKF executable
add_executable(letkf letkf.cpp)

# Link with required backend libraries for logging and configuration
target_link_libraries(letkf
    PRIVATE
        metada::traits
        metada::framework::adapters
        metada::framework::core
        metada::framework::interfaces
        metada::framework::runs
        metada::backends::common::utils
)

# Set output directory
set_target_properties(letkf
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Add format target to ensure consistent code style across the project
AddFormatTarget(letkf ${CMAKE_CURRENT_SOURCE_DIR})