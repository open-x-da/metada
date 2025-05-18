# Find required dependencies for MACOM
# Ensure FindMPI is robust. It should set MPI_FOUND and imported targets like MPI::MPI_C, MPI::MPI_Fortran.
metada_find_package(MPI REQUIRED) # Ensure this is effective
metada_find_package(HDF5 REQUIRED)
metada_find_package(netCDF REQUIRED)
metada_find_package(CURL REQUIRED)

# Create the MACOM backend library
add_library(macom_backends INTERFACE)
add_library(metada::backends::macom ALIAS macom_backends)

# set(MSYS2_PREFIX "/mingw64")
# set(ZLIB_LIB ${MSYS2_PREFIX}/lib/libz.dll.a)
# set(SZIP_LIB ${MSYS2_PREFIX}/lib/libsz.dll.a)
# set(MPI_LIB  ${MSYS2_PREFIX}/lib/libmsmpi.dll.a)

# set(COMMON_INCLUDE ${MSYS2_PREFIX}/include)

# CAT Preprocessing Flags
option(SEAICE_VARIABLE_SALINITY "Variable salinity" ON)
option(SEAICE_ITD "Sea Ice Thickness Distribution" ON)
option(SEAICE_ALLOW_FREEDRIFT "Enable free drift code" ON)
option(SEAICE_ALLOW_DYNAMICS "Allow sea ice dynamic code" ON)
option(SEAICE_ALLOW_MOM_ADVECTION "Enable advection terms" OFF)
option(SEAICE_ALLOW_BOTTOMDRAG "Use parameterization of grounding ice" ON)
option(SEAICE_CAP_SUBLIM "Used in thermodyn module" ON)
option(OPENACC "Enable OpenACC" OFF)
option(OPENACCGPU "Enable OpenACC GPU" OFF)
option(SeaiceDebug "Enable Seaice Debug" OFF)

set(DEF_FLAGS "")
if(SEAICE_VARIABLE_SALINITY)
    list(APPEND DEF_FLAGS "-DSEAICE_VARIABLE_SALINITY")
endif()
if(SEAICE_ITD)
    list(APPEND DEF_FLAGS "-DSEAICE_ITD")
endif()
if(SEAICE_ALLOW_FREEDRIFT)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_FREEDRIFT")
endif()
if(SEAICE_ALLOW_DYNAMICS)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_DYNAMICS")
endif()
if(SEAICE_ALLOW_MOM_ADVECTION)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_MOM_ADVECTION")
endif()
if(SEAICE_ALLOW_BOTTOMDRAG)
    list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_BOTTOMDRAG")
endif()
if(SEAICE_CAP_SUBLIM)
    list(APPEND DEF_FLAGS "-DSEAICE_CAP_SUBLIM")
endif()
if(OPENACC)
    list(APPEND DEF_FLAGS "-DOPENACC")
endif()
if(OPENACCGPU)
    list(APPEND DEF_FLAGS "-DOPENACCGPU")
endif()
if(SeaiceDebug)
    list(APPEND DEF_FLAGS "-DSeaiceDebug")
endif()

# set(CMAKE_Fortran_COMPILER mpifort)
# set(CMAKE_Fortran_FLAGS "-O3 -cpp ${DEF_FLAGS_STR}")


# Set include directories
target_include_directories(macom_backends 
    INTERFACE  
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${COMMON_INCLUDE}>
)

# Set headers as interface sources to show them in IDE
target_sources(macom_backends
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include/MACOMModel.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/MACOMState.hpp
)

# Link dependencies
target_link_libraries(macom_backends 
    INTERFACE
        HDF5
        # NetCDF
        # SZIP
        # ZLIB
        # CURL
        # MPI
)

# Add form# at target
AddFormatTarget(macom_backends ${CMAKE_CURRENT_SOURCE_DIR})


# Assume you have a list of sources for the macom backend library
set(MACOM_BACKEND_SOURCES
    src/MACOMFortranInterface.cpp
    # src/MACOMModel.cpp
    # ... other .cpp files for macom backend ...
)

# Assume you have a list of fortran sources
set(MACOM_FORTRAN_SOURCES
    fortran/macom_fortran_wrapper.f90
    # ... other .f90 files ...
)

# Create a library for the MACOM backend
add_library(macom_backend ${MACOM_BACKEND_SOURCES} ${MACOM_FORTRAN_SOURCES})

# Link macom_backend to MPI and other necessary libraries
target_link_libraries(macom_backend PUBLIC
    HDF5
    NetCDF
    CURL
    # MPI::MPI_C       # For C++ MPI calls if any, and for MACOMFortranInterface
    # MPI::MPI_Fortran # For the Fortran MPI calls in macom_fortran_wrapper.f90
    # ... other libraries macom_backend depends on ...
)

# Ensure this library is then linked by the forecast executable
# In applications/forecast/CMakeLists.txt:
# target_link_libraries(forecast PRIVATE macom_backend ...)

# Set target properties for macom_backend

# Include directories for macom_backend (both for its sources and for targets linking to it)
target_include_directories(macom_backend PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>  # For MACOMModel.hpp, etc.
    # Add other necessary include paths for MACOM's internal Fortran modules if they are
    # in different subdirectories and referenced by macom_fortran_wrapper.f90
    # e.g., $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/modules>
)

# Set Fortran compiler flags (including preprocessor definitions)
# CMAKE_Fortran_FLAGS will apply globally or to all Fortran targets if set here.
# It's often better to set per-target using target_compile_options or target_compile_definitions.
# However, if DEF_FLAGS are meant for all Fortran code in this subproject:
set_property(TARGET macom_backend PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/macom) # Output .mod files here

# Join DEF_FLAGS into a string for target_compile_options if needed for Fortran
string(REPLACE ";" " " DEF_FLAGS_STR "${DEF_FLAGS}")

# Apply compile options/definitions to the Fortran part of macom_backend
# This is a more modern way than setting CMAKE_Fortran_FLAGS directly for all targets.
# Note: target_compile_definitions is primarily for C/C++. For Fortran preprocessor, 
# passing them via compile options is common, or using a custom build step if complex.
# For simple -D flags, some compilers accept them via CMAKE_Fortran_FLAGS or target_compile_options.
# Example for gfortran if it takes -D flags directly in compile options:
target_compile_options(macom_backend PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${DEF_FLAGS_STR}>
)
# Alternatively, if your metada_find_package or toolchain handles Fortran flags:
# get_property(METADA_FORTRAN_FLAGS GLOBAL PROPERTY METADA_FORTRAN_FLAGS)
# target_compile_options(macom_backend PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${METADA_FORTRAN_FLAGS}>)

# Link macom_backend against its dependencies
# target_link_libraries(macom_backend PUBLIC
#     # MPI libraries for both C++ and Fortran
#     MPI::MPI_C       # For MACOMFortranInterface.cpp, MACOMModel.cpp
#     MPI::MPI_Fortran # For macom_fortran_wrapper.f90

#     # Other dependencies
#     # Using the imported targets from find_package is preferred
#     HDF5::HDF5       # Assuming HDF5 created HDF5::HDF5 (check your FindHDF5.cmake)
#                      # If it created HDF5::hdf5, use that. Or HDF5_LIBRARIES if not an imported target.
#     netCDF::netcdf_c # Or similar imported target from FindNetCDF.cmake
#     CURL::libcurl    # Or similar imported target from FindCURL.cmake

#     # If your FindXXX.cmake scripts only provide variables like XXX_LIBRARIES:
#     # ${HDF5_LIBRARIES}
#     # ${NETCDF_LIBRARIES}
#     # ${CURL_LIBRARIES}
# )

# Add format target for macom_backend sources
# Assuming AddFormatTarget is a custom function you have defined.
# It should probably take the source files list instead of CMAKE_CURRENT_SOURCE_DIR
# AddFormatTarget(macom_backend_sources ${MACOM_CPP_SOURCES} ${MACOM_FORTRAN_SOURCES})
# Or if it operates on a target:
AddFormatTarget(macom_backend ${CMAKE_CURRENT_SOURCE_DIR})


