# Find required dependencies for MACOM
# Ensure FindMPI is robust. It should set MPI_FOUND and imported targets like MPI::MPI_C, MPI::MPI_Fortran.
# MACOM mode selection - mutually exclusive options
option(MACOM_REGION_MODE "Build MACOM in region mode" OFF)
option(MACOM_GLOBAL_MODE "Build MACOM in global mode" OFF)

# Check if any MACOM mode is enabled
set(MACOM_MODE_ENABLED FALSE)
if(MACOM_REGION_MODE OR MACOM_GLOBAL_MODE)
    set(MACOM_MODE_ENABLED TRUE)
endif()

# Ensure only one mode is selected if any mode is enabled
if(MACOM_REGION_MODE AND MACOM_GLOBAL_MODE)
    message(FATAL_ERROR "MACOM_REGION_MODE and MACOM_GLOBAL_MODE are mutually exclusive. Please select only one.")
endif()

if(MACOM_MODE_ENABLED)
    metada_find_package(MPI REQUIRED)
    metada_find_package(HDF5 REQUIRED)
    metada_find_package(CURL REQUIRED)
else()
    # When no MACOM mode is enabled, make MPI optional
    metada_find_package(MPI)
    metada_find_package(HDF5)
    metada_find_package(CURL)
endif()

if(WIN32)
    # Find NetCDF Fortran for MSYS2/MINGW64
    set(NetCDF_Fortran_ROOT "C:/msys64/mingw64")
    set(NetCDF_Fortran_INCLUDE_DIRS "${NetCDF_Fortran_ROOT}/include")
    set(NetCDF_Fortran_LIBRARIES "${NetCDF_Fortran_ROOT}/lib/libnetcdff.dll.a")
else()
    # Linux/Unix platforms
    set(NetCDF_Fortran_INCLUDE_DIRS "/usr/include")
    set(NetCDF_Fortran_LIBRARIES "/usr/lib/x86_64-linux-gnu/libnetcdff.so")
    # Ensure variables are properly set
    if(NOT HDF5_LIBRARIES)
        set(HDF5_INCLUDE_DIRS /usr/include)
        set(HDF5_LIBRARIES /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so)
    endif()

    if(NOT CURL_LIBRARIES)
        set(CURL_INCLUDE_DIRS /usr/include)
        set(CURL_LIBRARIES /usr/lib/x86_64-linux-gnu/libcurl.so)
    endif()
endif()
metada_find_package(xtensor CONFIG)
metada_find_package(netcdf-cxx4)
metada_find_package(netCDF CONFIG)
metada_find_package(xtl CONFIG)
metada_find_package(nanoflann CONFIG)

# Create the MACOM backend interface library
add_library(macom_backends INTERFACE)
add_library(metada::backends::macom ALIAS macom_backends)

# Set include directories for the interface library
target_include_directories(macom_backends 
    INTERFACE  
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Set headers as interface sources to show them in IDE
target_sources(macom_backends
    INTERFACE
        MACOMModel.hpp
        MACOMGeometry.hpp
        MACOMGeometryIterator.hpp
        MACOMState.hpp
)

# Link dependencies for C++ interface # xtensor
target_link_libraries(macom_backends 
    INTERFACE
        metada::base
        NetCDF::CXX4
        nanoflann::nanoflann
)

# Only proceed with Fortran compilation if a mode is enabled
if(MACOM_MODE_ENABLED)
    # Set Fortran compiler
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

    # Set compiler flags
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -cpp -fno-tree-vectorize -fallow-argument-mismatch -Wno-missing-include-dirs")

    # Fortran compilation flags and macros
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

    # Auto-disable sea ice options in region mode
    if(MACOM_REGION_MODE)
        set(SEAICE_VARIABLE_SALINITY OFF)
        set(SEAICE_ITD OFF)
        set(SEAICE_ALLOW_FREEDRIFT OFF)
        set(SEAICE_ALLOW_DYNAMICS OFF)
        set(SEAICE_ALLOW_MOM_ADVECTION OFF)
        set(SEAICE_ALLOW_BOTTOMDRAG OFF)
        set(SEAICE_CAP_SUBLIM OFF)
        set(SeaiceDebug OFF)
        message(STATUS "Auto-disabled sea ice options in region mode")
    endif()

    # Set Fortran preprocessor macro definitions
    set(DEF_FLAGS "")

    # Sea ice options - available in both modes but can be disabled
    if(SEAICE_VARIABLE_SALINITY)
        list(APPEND DEF_FLAGS "-DSEAICE_VARIABLE_SALINITY")
        add_compile_definitions(SEAICE_VARIABLE_SALINITY)
    endif()
    if(SEAICE_ITD)
        list(APPEND DEF_FLAGS "-DSEAICE_ITD")
        add_compile_definitions(SEAICE_ITD)
    endif()
    if(SEAICE_ALLOW_FREEDRIFT)
        list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_FREEDRIFT")
        add_compile_definitions(SEAICE_ALLOW_FREEDRIFT)
    endif()
    if(SEAICE_ALLOW_DYNAMICS)
        list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_DYNAMICS")
        add_compile_definitions(SEAICE_ALLOW_DYNAMICS)
    endif()
    if(SEAICE_ALLOW_MOM_ADVECTION)
        list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_MOM_ADVECTION")
        add_compile_definitions(SEAICE_ALLOW_MOM_ADVECTION)
    endif()
    if(SEAICE_ALLOW_BOTTOMDRAG)
        list(APPEND DEF_FLAGS "-DSEAICE_ALLOW_BOTTOMDRAG")
        add_compile_definitions(SEAICE_ALLOW_BOTTOMDRAG)
    endif()
    if(SEAICE_CAP_SUBLIM)
        list(APPEND DEF_FLAGS "-DSEAICE_CAP_SUBLIM")
        add_compile_definitions(SEAICE_CAP_SUBLIM)
    endif()
    if(OPENACC)
        list(APPEND DEF_FLAGS "-DOPENACC")
        add_compile_definitions(OPENACC)
    endif()
    if(OPENACCGPU)
        list(APPEND DEF_FLAGS "-DOPENACCGPU")
        add_compile_definitions(OPENACCGPU)
    endif()
    if(SeaiceDebug)
        list(APPEND DEF_FLAGS "-DSeaiceDebug")
        add_compile_definitions(SeaiceDebug)
    endif()

    # Mode-specific flags
    if(MACOM_REGION_MODE)
        list(APPEND DEF_FLAGS "-DBARMOD")
        add_compile_definitions(BARMOD)
        message(STATUS "Building MACOM in REGION mode (sea ice auto-disabled)")
    else()
        message(STATUS "Building MACOM in GLOBAL mode")
    endif()

    string(REPLACE ";" " " DEF_FLAGS_STR "${DEF_FLAGS}")

    # Ensure Fortran module directory is explicitly set in CMakeLists.txt
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

    # Select source files based on mode
    if(MACOM_REGION_MODE)
        # Region mode: use src-region files
        file(GLOB SRC_FORTRAN_FILES
            "${CMAKE_CURRENT_SOURCE_DIR}/fortran/src-region/*.f90"
        )
        message(STATUS "Using region mode Fortran sources from src-region/")
    else()
        # Global mode: use src-globe files
        file(GLOB SRC_FORTRAN_FILES
            "${CMAKE_CURRENT_SOURCE_DIR}/fortran/src-globe/*.f90"
        )
        message(STATUS "Using global mode Fortran sources from src-globe/")
    endif()

    # Create Fortran implementation library, ensure macom_logger.f90 comes before macom_fortran_wrapper.f90
    add_library(macom_fortran STATIC
        fortran/macom_logger.f90           
        fortran/macom_fortran_wrapper.f90
        ${SRC_FORTRAN_FILES}
    )

    # Set Fortran module output directory
    set_property(TARGET macom_fortran PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

    # Add module directory to include paths
    target_include_directories(macom_fortran PUBLIC
        ${CMAKE_Fortran_MODULE_DIRECTORY}
        ${NetCDF_Fortran_INCLUDE_DIRS}
        $<IF:$<NOT:$<BOOL:${WIN32}>>,${HDF5_INCLUDE_DIRS};${CURL_INCLUDE_DIRS},>
    )

    # Set Fortran compiler flags and macros
    target_compile_definitions(macom_fortran PUBLIC
        $<$<COMPILE_LANGUAGE:Fortran>:${DEF_FLAGS}>
    )

    if(WIN32)
        target_link_libraries(macom_fortran PUBLIC HDF5 CURL)
    else()
        target_link_libraries(macom_fortran PUBLIC ${HDF5_LIBRARIES} ${CURL_LIBRARIES})
    endif()

    # Link Fortran dependencies
    target_link_libraries(macom_fortran PUBLIC
        ${NetCDF_Fortran_LIBRARIES}
        MPI::MPI_Fortran
    )

    # Create C++ implementation library that links with Fortran
    add_library(macom_impl STATIC
        src/MACOMFortranInterface.cpp
    )

    # Set include directories for implementation
    target_include_directories(macom_impl
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

    # Link C++ implementation with Fortran implementation and dependencies
    target_link_libraries(macom_impl PUBLIC
        macom_fortran
    )

    # Link interface with implementation
    target_link_libraries(macom_backends INTERFACE macom_impl)

else()
    # No MACOM mode enabled - only create C++ interface without Fortran
    message(STATUS "MACOM modes disabled - building C++ interface only (no Fortran compilation)")
    
    # Create C++ implementation library without Fortran
    add_library(macom_impl STATIC
        src/MACOMFortranInterface.cpp
    )

    # Set include directories for implementation
    target_include_directories(macom_impl
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

    # Link interface with implementation (no Fortran dependencies)
    target_link_libraries(macom_backends INTERFACE macom_impl)
endif()

# Add format targets
AddFormatTarget(macom_backends ${CMAKE_CURRENT_SOURCE_DIR})

# Fortran formatting configuration - Multiple options for different use cases:

# Option 1: Format all Fortran files (SLOW - not recommended)
# AddFortranFormatTarget(macom_backends ${CMAKE_CURRENT_SOURCE_DIR})

# Option 2: Exclude large directories from formatting (MEDIUM SPEED)
# AddFortranFormatTarget(macom_backends ${CMAKE_CURRENT_SOURCE_DIR} EXCLUDE_DIRS "fortran/src-region;fortran/src-globe")

# Option 3: Format only specific files (FASTEST - RECOMMENDED)
# Using relative paths from CMAKE_CURRENT_SOURCE_DIR
AddFortranFormatTargetSelective(macom_backends FILES "fortran/macom_logger.f90;fortran/macom_fortran_wrapper.f90")

# Alternative: Using absolute paths (more explicit)
# AddFortranFormatTargetSelective(macom_backends FILES "${CMAKE_CURRENT_SOURCE_DIR}/fortran/macom_logger.f90;${CMAKE_CURRENT_SOURCE_DIR}/fortran/macom_fortran_wrapper.f90")

# Option 4: Format specific directories only
# AddFortranFormatTargetSelective(macom_backends DIRS "fortran/wrapper;fortran/utils")

# Option 5: Disable Fortran formatting entirely
# (No AddFortranFormatTarget call)



