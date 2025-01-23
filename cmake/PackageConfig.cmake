# Package configuration lists
set(METADA_REQUIRED_PACKAGES
    CACHE INTERNAL "List of required packages"
    FORCE
)

set(METADA_OPTIONAL_PACKAGES
    CACHE INTERNAL "List of optional packages"
    FORCE
)

# Function to register required package
function(register_required_package package)
    list(APPEND METADA_REQUIRED_PACKAGES ${package})
    set(METADA_REQUIRED_PACKAGES ${METADA_REQUIRED_PACKAGES}
        CACHE INTERNAL "List of required packages"
        FORCE
    )
endfunction()

# Function to register optional package
function(register_optional_package package)
    list(APPEND METADA_OPTIONAL_PACKAGES ${package})
    set(METADA_OPTIONAL_PACKAGES ${METADA_OPTIONAL_PACKAGES}
        CACHE INTERNAL "List of optional packages"
        FORCE
    )
endfunction() 