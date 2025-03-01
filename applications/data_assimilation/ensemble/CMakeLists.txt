# Add LETKF executable
add_executable(letkf letkf.cpp)

# Link with required backend libraries for logging and configuration
target_link_libraries(letkf
    PRIVATE
        metada::framework::common
        metada::framework::runs
        metada::framework::traits
        metada::backends::common::utils::logger
        metada::backends::common::utils::config
)

# Add format target to ensure consistent code style across the project
AddFormatTarget(letkf ${CMAKE_CURRENT_SOURCE_DIR})