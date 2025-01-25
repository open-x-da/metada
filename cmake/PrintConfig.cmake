include(print/Headers)
include(print/PackageInfo)
include(print/CompilerInfo)
include(print/BuildInfo)
include(print/SystemInfo)

function(print_configuration_summary)
    print_section_header("Build Configuration")
    print_compiler_info()
    print_build_info()
    print_system_info()
    
    print_section_header("Required Packages")
    foreach(pkg IN LISTS METADA_REQUIRED_PACKAGES)
        print_package_info(${pkg})
    endforeach()
    
    print_section_header("Optional Packages")
    foreach(pkg IN LISTS METADA_OPTIONAL_PACKAGES)
        print_package_info(${pkg})
    endforeach()
endfunction()

print_configuration_summary()