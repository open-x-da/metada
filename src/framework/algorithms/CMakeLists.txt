# Create algorithms interface library
add_library(algorithms INTERFACE)

add_library(metada::framework::algorithms ALIAS algorithms)

# Set include directories
target_include_directories(algorithms
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link dependencies
target_link_libraries(algorithms
    INTERFACE metada::framework::adapters
    INTERFACE metada::base
    INTERFACE metada::traits
    INTERFACE Eigen3::Eigen
)

# Add format target for all algorithm headers
AddFormatTarget(algorithms ${CMAKE_CURRENT_SOURCE_DIR}) 