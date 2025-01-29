# Add main executable
add_executable(metada main.cpp)

# Link required libraries
target_link_libraries(metada
    PRIVATE
        glog::glog
)

# CUDA configuration
if(USE_CUDA)
    target_link_libraries(metada
        PRIVATE
            cuda_settings
    )
    target_include_directories(metada
        PRIVATE
            ${CUDAToolkit_INCLUDE_DIRS}
    )
else()
    target_compile_definitions(metada
        PRIVATE
            NO_CUDA
    )
endif()

# Set C++ standard
set_target_properties(metada
    PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)

# Clean coverage data before build
CleanCoverage(metada)

# Add format target
add_format_target(metada ${CMAKE_CURRENT_SOURCE_DIR})

# Install rules
install(TARGETS metada
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
) 