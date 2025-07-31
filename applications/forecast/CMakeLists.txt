# Add MPI support option
option(USE_MPI "Enable MPI support for forecast application" OFF)

if(USE_MPI)
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        add_definitions(-DUSE_MPI)
        include_directories(${MPI_INCLUDE_PATH})
        message(STATUS "MPI support enabled for forecast application")
        message(STATUS "MPI_CXX_INCLUDE_DIRS: ${MPI_CXX_INCLUDE_DIRS}")
        message(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
    else()
        message(FATAL_ERROR "MPI support requested but MPI not found")
    endif()
else()
    message(STATUS "MPI support disabled for forecast application - building serial version")
endif()

# Create forecast executable
add_executable(forecast forecast.cpp)

# Link with all required libraries
target_link_libraries(forecast
    PRIVATE
        metada::base
        metada::framework::adapters
        metada::framework::runs
        metada::traits
        metada::backends::common::utils
        metada::backends::common::observation
        metada::backends::simple
        # metada::backends::wrf
        metada::backends::macom
        $<$<BOOL:${USE_MPI}>:${MPI_CXX_LIBRARIES}>
)

# Set output directory
set_target_properties(forecast
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Add format target to ensure consistent code style
AddFormatTarget(forecast ${CMAKE_CURRENT_SOURCE_DIR})