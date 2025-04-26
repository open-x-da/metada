# Create forecast executable
add_executable(forecast forecast.cpp)

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

# Add format target to ensure consistent code style
AddFormatTarget(forecast ${CMAKE_CURRENT_SOURCE_DIR})