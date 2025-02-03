# Add LETKF executable
add_executable(letkf letkf.cc)

# Link only with framework logger
target_link_libraries(letkf
    PRIVATE
        metada::framework::tools::logger
)

# Add format target to ensure consistent code style across the project
AddFormatTarget(letkf ${CMAKE_CURRENT_SOURCE_DIR})