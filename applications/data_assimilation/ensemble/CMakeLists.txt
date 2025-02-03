# Add LETKF executable
add_executable(letkf letkf.cc)

# Link required libraries
target_link_libraries(letkf
    PRIVATE
        metada::framework::tools::logger
)

if(USE_GLOG)
    target_link_libraries(letkf
        PRIVATE
            metada::backends::tools::logger
    )
endif()

# Add format target to ensure consistent code style across the project
AddFormatTarget(letkf ${CMAKE_CURRENT_SOURCE_DIR})