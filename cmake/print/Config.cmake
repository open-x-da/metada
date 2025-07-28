include(print/Formatting)
include(print/PackageInfo)
include(print/CompilerInfo)
include(print/BuildInfo)
include(print/SystemInfo)

function(metada_project_summary)
    print_header("METADA Configuration Summary")
    
    # System information
    print_subheader("System Information")
    print_system_info()

    # Compiler information
    print_subheader("Compiler Setup")
    print_compiler_info()
    
    # Build information
    print_subheader("Build Configuration") 
    print_build_info()

    # Required packages
    print_subheader("Required Dependencies")
    foreach(pkg IN LISTS METADA_REQUIRED_PACKAGES)
        print_package_info(${pkg})
    endforeach()
    
    # Optional packages
    print_subheader("Optional Dependencies")
    foreach(pkg IN LISTS METADA_OPTIONAL_PACKAGES)
        print_package_info(${pkg})
    endforeach()
    
    message(STATUS "\n") # Add final newline
endfunction()