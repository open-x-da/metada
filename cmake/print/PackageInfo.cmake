# Print package information
include(print/Headers)
include(package/Path)

function(print_package_info package)
    if(NOT ${package}_FOUND)
        message(STATUS "  - ${package}: Not found")
        return()
    endif()

    # Print basic info
    if(${package}_VERSION)
        message(STATUS "  - ${package}: Found (${${package}_VERSION})")
    else()
        message(STATUS "  - ${package}: Found")
    endif()

    # Print path if available
    get_package_path(${package} path)
    if(path)
        message(STATUS "    Path: ${path}")
    endif()

    # Print components
    get_package_components(${package} components)
    if(components)
        message(STATUS "    Components:")
        foreach(component ${components})
            if(${package}_${component}_FOUND)
                message(STATUS "      - ${component}: Found")
            else()
                message(STATUS "      - ${component}: Not found")
            endif()
        endforeach()
    endif()
endfunction() 