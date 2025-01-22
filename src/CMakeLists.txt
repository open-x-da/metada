# Add main executable
add_executable(metada main.cpp)

# Link required libraries
target_link_libraries(metada 
    PRIVATE 
        CUDA::cudart
        glog::glog
)

# Set include directories
target_include_directories(metada
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Set C++ standard
set_target_properties(metada PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

# Install rules
install(TARGETS metada
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
) 