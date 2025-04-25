# Add the WRF Backend sources
set(WRF_BACKEND_SOURCES
    Geometry.cpp
    GeometryIterator.cpp
    State.cpp
    Model.cpp
)

# Find required dependencies
find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(netcdf-cxx4 REQUIRED)

# Create the WRF backend library
add_library(metada_wrf_backend ${WRF_BACKEND_SOURCES})
add_library(metada::backends::wrf ALIAS metada_wrf_backend)

# Set include directories
target_include_directories(metada_wrf_backend 
    PUBLIC 
        ${CMAKE_SOURCE_DIR}/src
    SYSTEM PUBLIC
        ${xtl_INCLUDE_DIRS}
        ${xtensor_INCLUDE_DIRS}
        ${netcdf-cxx4_INCLUDE_DIRS}
)

# Link dependencies
target_link_libraries(metada_wrf_backend 
    PUBLIC
        ${netcdf-cxx4_LIBRARIES}
        xtl
        xtensor
        metada::framework::core
)

# Add format target
AddFormatTarget(metada_wrf_backend ${CMAKE_CURRENT_SOURCE_DIR}) 