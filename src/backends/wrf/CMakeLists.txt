# Add the WRF Backend sources
set(WRF_BACKEND_SOURCES
    WRFGeometry.cpp
    WRFState.cpp
    WRFModel.cpp
)

# Find required dependencies
metada_find_package(xtl CONFIG)
metada_find_package(xtensor CONFIG)
metada_find_package(netcdf-cxx4)

# Create the WRF backend library
add_library(wrf_backends ${WRF_BACKEND_SOURCES})
add_library(metada::backends::wrf ALIAS wrf_backends)

# Set include directories
target_include_directories(wrf_backends 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Link dependencies
target_link_libraries(wrf_backends 
    PUBLIC
        NetCDF::CXX4
        xtl
        xtensor
        metada::framework::core
)

# Add format target
AddFormatTarget(wrf_backends ${CMAKE_CURRENT_SOURCE_DIR}) 
