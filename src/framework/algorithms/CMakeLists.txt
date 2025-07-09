# Create algorithms interface library
add_library(algorithms INTERFACE)

add_library(metada::framework::algorithms ALIAS algorithms)

# Set include directories
target_include_directories(algorithms
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
)

# List sources for IDE/documentation
target_sources(algorithms
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/ETKF.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/LETKF.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Metrics.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ProgressBar.hpp
)

# Link dependencies
target_link_libraries(algorithms
    INTERFACE metada::framework::adapters
    INTERFACE metada::base
    INTERFACE metada::traits
    INTERFACE Eigen3::Eigen
)

AddFormatTarget(algorithms ${CMAKE_CURRENT_SOURCE_DIR}) 