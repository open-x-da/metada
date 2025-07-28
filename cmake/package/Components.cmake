# Handle package components tracking
function(register_package_components package)
    cmake_parse_arguments(ARG "" "" "COMPONENTS" ${ARGN})
    if(ARG_COMPONENTS)
        set_property(GLOBAL APPEND PROPERTY ${package}_COMPONENTS ${ARG_COMPONENTS})
    endif()
endfunction()

function(get_package_components package out_var)
    get_property(components GLOBAL PROPERTY ${package}_COMPONENTS)
    set(${out_var} "${components}" PARENT_SCOPE)
endfunction() 