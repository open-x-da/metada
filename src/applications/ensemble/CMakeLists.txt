# Add LETKF executable
add_executable(letkf letkf.cc)

# Link required libraries
target_link_libraries(letkf
    PRIVATE
        metada_core
        glog::glog
)

# Set C++ standard
set_target_properties(letkf
    PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)

# Add format target
add_format_target(letkf ${CMAKE_CURRENT_SOURCE_DIR})

# Install rules
install(TARGETS letkf
    RUNTIME DESTINATION bin
) 