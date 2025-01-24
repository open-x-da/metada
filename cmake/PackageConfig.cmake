# Package configuration lists
set(METADA_REQUIRED_PACKAGES
    CACHE INTERNAL "List of required packages"
    FORCE
)

set(METADA_OPTIONAL_PACKAGES
    CACHE INTERNAL "List of optional packages"
    FORCE
)

# Function to find and register a package
function(metada_find_package package)
    cmake_parse_arguments(ARG
        "OPTIONAL;QUIET"
        "CONDITION"
        "COMPONENTS"
        ${ARGN}
    )

    # Check condition if specified
    if(DEFINED ARG_CONDITION AND NOT ${ARG_CONDITION})
        return()
    endif()

    # Prepare find_package arguments
    set(find_args ${package})
    if(ARG_COMPONENTS)
        list(APPEND find_args COMPONENTS ${ARG_COMPONENTS})
    endif()
    if(ARG_QUIET)
        list(APPEND find_args QUIET)
    endif()
    if(NOT ARG_OPTIONAL)
        list(APPEND find_args REQUIRED)
    endif()

    # Find package
    find_package(${find_args})

    # Set _FOUND variable in cache
    set(${package}_FOUND ${${package}_FOUND} 
        CACHE INTERNAL "${package} found status"
        FORCE
    )
    
    # Set version in cache if it exists
    if(DEFINED ${package}_VERSION)
        set(${package}_VERSION "${${package}_VERSION}"
            CACHE INTERNAL "${package} version"
            FORCE
        )
    endif()

    # Register package
    if(ARG_OPTIONAL)
        if(NOT "${package}" IN_LIST METADA_OPTIONAL_PACKAGES)
            list(APPEND METADA_OPTIONAL_PACKAGES ${package})
            set(METADA_OPTIONAL_PACKAGES ${METADA_OPTIONAL_PACKAGES}
                CACHE INTERNAL "List of optional packages"
                FORCE
            )
        endif()
    else()
        if(NOT "${package}" IN_LIST METADA_REQUIRED_PACKAGES)
            list(APPEND METADA_REQUIRED_PACKAGES ${package})
            set(METADA_REQUIRED_PACKAGES ${METADA_REQUIRED_PACKAGES}
                CACHE INTERNAL "List of required packages"
                FORCE
            )
        endif()
    endif()
endfunction() 