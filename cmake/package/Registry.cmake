# Package registry management
set(METADA_REQUIRED_PACKAGES "" CACHE INTERNAL "List of required packages")
set(METADA_OPTIONAL_PACKAGES "" CACHE INTERNAL "List of optional packages")

function(register_package package)
    cmake_parse_arguments(ARG "OPTIONAL" "" "" ${ARGN})
    
    if(ARG_OPTIONAL)
        list(APPEND METADA_OPTIONAL_PACKAGES ${package})
        set(METADA_OPTIONAL_PACKAGES ${METADA_OPTIONAL_PACKAGES} CACHE INTERNAL "" FORCE)
    else()
        list(APPEND METADA_REQUIRED_PACKAGES ${package})
        set(METADA_REQUIRED_PACKAGES ${METADA_REQUIRED_PACKAGES} CACHE INTERNAL "" FORCE)
    endif()
endfunction() 