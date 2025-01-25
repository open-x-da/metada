include(package/Components)
include(package/Registry)
include(package/Status)

# Function to find and register a package
function(metada_find_package package)
    cmake_parse_arguments(ARG
        "OPTIONAL;QUIET"
        "CONDITION"
        "COMPONENTS"
        ${ARGN}
    )

    # Early return if condition not met
    if(DEFINED ARG_CONDITION AND NOT ${ARG_CONDITION})
        return()
    endif()

    # Find package
    set(find_args ${package})
    if(ARG_COMPONENTS)
        list(APPEND find_args COMPONENTS ${ARG_COMPONENTS})
    endif()
    if(NOT ARG_OPTIONAL)
        list(APPEND find_args REQUIRED)
    endif()
    if(ARG_QUIET)
        list(APPEND find_args QUIET)
    endif()
    find_package(${find_args})

    # Register package info
    register_package_components(${package} COMPONENTS ${ARG_COMPONENTS})
    if(ARG_OPTIONAL)
        register_package(${package} OPTIONAL)
    else()
        register_package(${package})
    endif()
    set_package_status(${package} COMPONENTS ${ARG_COMPONENTS})
endfunction() 