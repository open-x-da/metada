# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(obs_processing obs_processing.cpp)

# Link with required libraries
target_link_libraries(obs_processing
    PRIVATE
        metada::base
        metada::framework::adapters
        metada::framework::runs
        metada::traits
        metada::backends::common::utils
        metada::backends::common::io
)

# Set output directory
set_target_properties(obs_processing
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Add format target to ensure consistent code style
AddFormatTarget(obs_processing ${CMAKE_CURRENT_SOURCE_DIR})