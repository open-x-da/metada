# Helper function to print a section header
function(_print_section_header text)
    message(STATUS "")
    message(STATUS "${text}")
    string(LENGTH "${text}" text_length)
    string(REPEAT "-" ${text_length} separator)
    message(STATUS "${separator}")
    message(STATUS "")
endfunction()

# Helper function to print package info
function(_print_package_info package_name)
    if(${package_name}_FOUND)
        set(version "${${package_name}_VERSION}")
        if(version)
            message(STATUS "  - ${package_name}: Found (${version})")
        else()
            message(STATUS "  - ${package_name}: Found")
        endif()
    else()
        message(STATUS "  - ${package_name}: Not found")
    endif()
endfunction()

# Main configuration printing function
function(print_configuration_summary)
    _print_section_header("Build Configuration")
    
    # Print compiler info
    message(STATUS "Compiler Information:")
    message(STATUS "  - C compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
    message(STATUS "  - C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "  - Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
    
    # Print build type
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
    
    # Print required packages
    _print_section_header("Required Packages")
    foreach(pkg IN LISTS METADA_REQUIRED_PACKAGES)
        _print_package_info(${pkg})
    endforeach()
    
    # Print optional packages
    _print_section_header("Optional Packages")
    foreach(pkg IN LISTS METADA_OPTIONAL_PACKAGES)
        _print_package_info(${pkg})
    endforeach()
endfunction() 