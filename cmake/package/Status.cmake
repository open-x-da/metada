# Package status handling
function(set_package_status package)
    cmake_parse_arguments(ARG "" "VERSION" "COMPONENTS" ${ARGN})
    
    string(TOUPPER ${package} PACKAGE_UPPER)
    
    # Handle package found status
    if(DEFINED ${PACKAGE_UPPER}_FOUND)
        set(${package}_FOUND ${${PACKAGE_UPPER}_FOUND} CACHE INTERNAL "" FORCE)
    elseif(DEFINED ${package}_FOUND)
        set(${package}_FOUND ${${package}_FOUND} CACHE INTERNAL "" FORCE)
    endif()
    
    # Handle version
    if(DEFINED ${PACKAGE_UPPER}_VERSION)
        set(${package}_VERSION "${${PACKAGE_UPPER}_VERSION}" CACHE INTERNAL "" FORCE)
    elseif(DEFINED ${package}_VERSION)
        set(${package}_VERSION "${${package}_VERSION}" CACHE INTERNAL "" FORCE)
    endif()
    
    # Handle components
    foreach(component ${ARG_COMPONENTS})
        string(TOUPPER ${component} COMPONENT_UPPER)
        if(DEFINED ${PACKAGE_UPPER}_${COMPONENT_UPPER}_FOUND)
            set(${package}_${component}_FOUND ${${PACKAGE_UPPER}_${COMPONENT_UPPER}_FOUND}
                CACHE INTERNAL "" FORCE)
        elseif(DEFINED ${package}_${component}_FOUND)
            set(${package}_${component}_FOUND ${${package}_${component}_FOUND}
                CACHE INTERNAL "" FORCE)
        endif()
    endforeach()
endfunction() 