# Add LETKF executable
add_executable(letkf letkf.cpp)

# Link only with backends logger
target_link_libraries(letkf
    PRIVATE
        metada::backends::tools::logger
        metada::backends::tools::config
)

# Add format target to ensure consistent code style across the project
AddFormatTarget(letkf ${CMAKE_CURRENT_SOURCE_DIR})