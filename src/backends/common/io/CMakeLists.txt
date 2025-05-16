metada_find_package(bufr CONFIG)

# Create the common io library (both interface and implementation)
add_library(common_io STATIC)
add_library(metada::backends::common::io ALIAS common_io)

# Add Fortran source files
target_sources(common_io
    PRIVATE
        readpb.f90
        virtmp.f90
        fortran_bufr_wrapper.f90
)

# Set Fortran fixed form source format
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

# Add C++ header files for IDE integration
set(HEADER_FILES
    BufrFortranAPI.h
    BufrObsIO.hpp
)

# Add include directories for headers
target_include_directories(common_io
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Make sure readpb.prm is found by the Fortran compiler
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/readpb.prm
    ${CMAKE_CURRENT_BINARY_DIR}/readpb.prm
    COPYONLY
)

# Link with required libraries
target_link_libraries(common_io
    PUBLIC
        bufr::bufr_4
        metada::base
)

# Add format target
AddFormatTarget(common_io ${CMAKE_CURRENT_SOURCE_DIR})
