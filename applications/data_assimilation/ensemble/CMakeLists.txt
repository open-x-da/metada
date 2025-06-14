# Add ETKF executable
add_executable(etkf main.cpp)

# Link with required backend libraries for logging and configuration
target_link_libraries(etkf
    PRIVATE
        metada::traits
        metada::framework::adapters
        metada::base
        metada::framework::runs
        metada::backends::common::utils
        metada::framework::algorithms
        metada::backends::gmock
)

# Set output directory
set_target_properties(etkf
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Add format target to ensure consistent code style across the project
AddFormatTarget(etkf ${CMAKE_CURRENT_SOURCE_DIR})