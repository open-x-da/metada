# Find BUFR library (optional)
set(BUFR_AVAILABLE FALSE)
if(WITH_BUFR)
    metada_find_package(bufr OPTIONAL CONFIG)
    # Check if BUFR target exists (more reliable than _FOUND variable)
    if(TARGET bufr::bufr_4)
        set(BUFR_AVAILABLE TRUE)
    else()
        message(WARNING "WITH_BUFR is enabled but bufr package not found. BUFR functionality will be disabled.")
        set(WITH_BUFR FALSE)
    endif()
endif()

# Enable Fortran language if BUFR is available
if(WITH_BUFR)
    enable_language(Fortran)
endif()

# Create the common io library (both interface and implementation)
add_library(common_io STATIC)
add_library(metada::backends::common::io ALIAS common_io)

# Always add stub source file (required for STATIC library)
target_sources(common_io
    PRIVATE
        stub.cpp
)

# Add Fortran source files (only when BUFR is enabled)
if(BUFR_AVAILABLE)
    target_sources(common_io
        PRIVATE
            readpb.f90
            virtmp.f90
            fortran_bufr_wrapper.f90
    )
endif()

# Set Fortran fixed form source format (only when BUFR is enabled)
if(BUFR_AVAILABLE)
    set_source_files_properties(
        readpb.f90
        virtmp.f90
        PROPERTIES
            Fortran_FORMAT FIXED
    )

    # Add compiler-specific flags for readpb.f90 to address stack size warnings
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set_source_files_properties(readpb.f90 PROPERTIES COMPILE_FLAGS "-frecursive")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set_source_files_properties(readpb.f90 PROPERTIES COMPILE_FLAGS "-recursive")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang|Clang")
        set_source_files_properties(readpb.f90 PROPERTIES COMPILE_FLAGS "-frecursive")
    endif()

    # Set Fortran free form source format
    set_source_files_properties(
        fortran_bufr_wrapper.f90
        PROPERTIES
            Fortran_FORMAT FREE
    )
endif()

# Add C++ header files for IDE integration
set(HEADER_FILES
    BufrFortranAPI.h
    BufrObsIO.hpp
)

# Add include directories for headers
target_include_directories(common_io
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/backends/common/obsoperator>
)

# Make sure readpb.prm is found by the Fortran compiler (only when BUFR is enabled)
if(BUFR_AVAILABLE)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/readpb.prm
        ${CMAKE_CURRENT_BINARY_DIR}/readpb.prm
        COPYONLY
    )
endif()

# Link with required libraries
target_link_libraries(common_io
    PUBLIC
        metada::base
)

# Link BUFR library only when enabled and found
if(BUFR_AVAILABLE)
    target_link_libraries(common_io
        PUBLIC
            bufr::bufr_4
    )
endif()

# Add format target
AddFormatTarget(common_io ${CMAKE_CURRENT_SOURCE_DIR})
