# Create algorithms interface library
add_library(algorithms STATIC
    LETKF.cpp
)

add_library(metada::framework::algorithms ALIAS algorithms)

# Set include directories
target_include_directories(algorithms
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

# List sources for IDE/documentation (optional, for header-only)
target_sources(algorithms
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/LETKF.hpp
)

# Link dependencies
target_link_libraries(algorithms
    PUBLIC metada::framework::adapters
    PUBLIC metada::base
    PUBLIC metada::traits
    PUBLIC Eigen3::Eigen
)

AddFormatTarget(algorithms ${CMAKE_CURRENT_SOURCE_DIR}) 