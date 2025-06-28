# Add ETKF executable
add_executable(etkf etkf.cpp)

# Add LETKF executable
add_executable(letkf letkf.cpp)

# Link with required backend libraries for logging and configuration
target_link_libraries(etkf
    PRIVATE
        metada::traits
        metada::framework::adapters
        metada::base
        metada::framework::runs
        metada::backends::common::utils
        metada::backends::common::observation
        metada::framework::algorithms
        metada::backends::wrf
)

target_link_libraries(letkf
    PRIVATE
        metada::traits
        metada::framework::adapters
        metada::base
        metada::framework::runs
        metada::backends::common::utils
        metada::backends::common::observation
        metada::framework::algorithms
        metada::backends::wrf
)

# Set output directory
set_target_properties(etkf
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(letkf
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Add format target to ensure consistent code style across the project
AddFormatTarget(etkf ${CMAKE_CURRENT_SOURCE_DIR})
AddFormatTarget(letkf ${CMAKE_CURRENT_SOURCE_DIR})