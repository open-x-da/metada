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

    # Find package with its components
    find_package(${find_args})

    # Convert package name to uppercase for compatibility with some CMake variables
    string(TOUPPER ${package} PACKAGE_UPPER)
    
    # Set _FOUND variable in cache if it exists
    if(DEFINED ${PACKAGE_UPPER}_FOUND)
        set(${package}_FOUND ${${PACKAGE_UPPER}_FOUND} 
            CACHE INTERNAL "${package} found status"
            FORCE
        )
    elseif(DEFINED ${package}_FOUND)
        set(${package}_FOUND ${${package}_FOUND} 
            CACHE INTERNAL "${package} found status"
            FORCE
        )
    endif()
    
    # Set version in cache if it exists
    if(DEFINED ${PACKAGE_UPPER}_VERSION)
        set(${package}_VERSION "${${PACKAGE_UPPER}_VERSION}"
            CACHE INTERNAL "${package} version"
            FORCE
        )
    elseif(DEFINED ${package}_VERSION)
        set(${package}_VERSION "${${package}_VERSION}"
            CACHE INTERNAL "${package} version"
            FORCE
        )
    endif()

    # Set component found status in cache
    if(ARG_COMPONENTS)
        foreach(component ${ARG_COMPONENTS})
            string(TOUPPER ${component} COMPONENT_UPPER)
            # Check both uppercase and regular variable names
            if(DEFINED ${PACKAGE_UPPER}_${COMPONENT_UPPER}_FOUND)
                set(${package}_${component}_FOUND ${${PACKAGE_UPPER}_${COMPONENT_UPPER}_FOUND}
                    CACHE INTERNAL "${package} ${component} component found status"
                    FORCE
                )
            elseif(DEFINED ${package}_${component}_FOUND)
                set(${package}_${component}_FOUND ${${package}_${component}_FOUND}
                    CACHE INTERNAL "${package} ${component} component found status"
                    FORCE
                )
            endif()
        endforeach()
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