# Find required dependencies
metada_find_package(xtl CONFIG)
metada_find_package(xtensor CONFIG)
metada_find_package(netcdf-cxx4)

# Create the WRF backend library
add_library(wrf_backends INTERFACE)
add_library(metada::backends::wrf ALIAS wrf_backends)

# Set include directories
target_include_directories(wrf_backends 
    INTERFACE 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Link dependencies
target_link_libraries(wrf_backends 
    INTERFACE
        NetCDF::CXX4
        xtensor
)

# Add format target
AddFormatTarget(wrf_backends ${CMAKE_CURRENT_SOURCE_DIR}) 
