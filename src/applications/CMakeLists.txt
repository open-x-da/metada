# List all application sources
set(APPLICATION_SOURCES
    forecast.cpp
)

# Create forecast executable
add_executable(forecast ${APPLICATION_SOURCES})

# Set include directories
target_include_directories(forecast
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
)

# Link with required libraries
target_link_libraries(forecast
    PRIVATE
        metada::framework::core
        metada::framework::adapters
        metada::framework::runs
        metada::traits
        metada::backends::common::utils
        metada::backends::wrf
)

# Set output directory
set_target_properties(forecast
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Add format target
AddFormatTarget(forecast ${CMAKE_CURRENT_SOURCE_DIR}) 