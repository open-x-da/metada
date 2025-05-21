# Helper functions for setting consistent target properties

# Function to set C++20 standard on a target
function(metada_set_cpp20_for_target target)
    # Set C++20 standard
    set_target_properties(${target} PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
    )
    
    # Add compiler feature requirements
    target_compile_features(${target} PUBLIC cxx_std_20)
    
    # Add specific compiler flags for C++20 concepts if needed
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
            target_compile_options(${target} PRIVATE -fconcepts)
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10.0)
            target_compile_options(${target} PRIVATE -fconcepts-ts)
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        target_compile_options(${target} PRIVATE /Zc:__cplusplus /permissive-)
    endif()
endfunction()

# Function to set common properties for all library targets
function(metada_setup_library_target target)
    # Set C++20 standard
    metada_set_cpp20_for_target(${target})
    
    # Set output directories
    set_target_properties(${target} PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    
    # Include what you use - ensures clean headers
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set_target_properties(${target} PROPERTIES
            CXX_INCLUDE_WHAT_YOU_USE "iwyu"
        )
    endif()
    
    # Add other common properties as needed
    
    # Print setup confirmation
    message(STATUS "Set up C++20 properties for target: ${target}")
endfunction()

# Function to set common properties for executable targets
function(metada_setup_executable_target target)
    # Set C++20 standard
    metada_set_cpp20_for_target(${target})
    
    # Set output directories
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    
    # Print setup confirmation
    message(STATUS "Set up C++20 properties for executable: ${target}")
endfunction() 