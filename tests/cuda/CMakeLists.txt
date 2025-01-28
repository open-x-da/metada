if(ENABLE_CUDA_ACCELERATION)
    # Device Query Test
    add_executable(device_query_test device_query_test.cpp)
    
    # Set CUDA C++11 standard for device_query_test
    # This is required because:
    # 1. CUDA C++11 provides modern features needed for device query implementation
    # 2. Ensures compatibility with CUDA toolkit features and APIs
    # 3. Matches the minimum C++ standard used in the rest of the project
    target_compile_features(device_query_test PRIVATE cuda_std_11)

    # Configure CUDA-specific properties
    # CUDA_SEPARABLE_COMPILATION enables separate compilation of CUDA code
    # which improves build times and allows for better optimization
    # LANGUAGE CUDA tells CMake this is a CUDA project requiring NVCC compiler
    set_target_properties(device_query_test PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        LANGUAGE CUDA
    )

    target_link_libraries(device_query_test 
        PRIVATE 
            cuda_settings
    )

    target_include_directories(device_query_test
        PRIVATE
            ${CMAKE_SOURCE_DIR}/src
    )

    AddTests(device_query_test)

    # Add format target
    add_format_target(device_query_test ${CMAKE_CURRENT_SOURCE_DIR})

    # CUDA Test Suite
    add_executable(cuda_test cuda_test.cpp)
    target_link_libraries(cuda_test
        PRIVATE
            cuda_settings
            glog::glog
    )

    AddTests(cuda_test)

    # Add format target
    add_format_target(cuda_test ${CMAKE_CURRENT_SOURCE_DIR})
endif() 