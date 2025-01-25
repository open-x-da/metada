# Print package information
include(package/Path)

function(print_package_info package)
    if(NOT ${package}_FOUND)
        print_item("${package}: Not found")
        return()
    endif()

    # Print basic info
    if(${package}_VERSION)
        print_item("${package}: Found (version ${${package}_VERSION})")
    else()
        print_item("${package}: Found")
    endif()

    # Print path if available
    get_package_path(${package} path)
    if(path)
        print_subitem("Path: ${path}")
    endif()

    # Print components
    get_package_components(${package} components)
    if(components)
        print_subitem("Components:")
        foreach(component ${components})
            if(${package}_${component}_FOUND)
                message(STATUS "      - ${component}: Found")
            else()
                message(STATUS "      - ${component}: Not found")
            endif()
        endforeach()
    endif()
endfunction() 