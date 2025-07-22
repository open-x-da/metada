# Add ETKF executable
add_executable(etkf etkf.cpp)

# Add LETKF executable
add_executable(letkf letkf.cpp)

# Add Particle Filter executable
add_executable(particle_filter particle_filter.cpp)

# Add EnKF executable
add_executable(enkf enkf.cpp)

# Add LWEnKF executable
add_executable(lwenkf lwenkf.cpp)

# Add A4DEnVar executable
add_executable(a4denvar a4denvar.cpp)

# Link with required backend libraries for logging and configuration
target_link_libraries(etkf
    PRIVATE
        metada::traits
        metada::framework::adapters
        metada::base
        metada::framework::runs
        metada::backends::common::utils
        metada::backends::common::observation
        metada::framework::algorithms
        metada::backends::simple
)

target_link_libraries(letkf
    PRIVATE
        metada::traits
        metada::framework::adapters
        metada::base
        metada::framework::runs
        metada::backends::common::utils
        metada::backends::common::observation
        metada::framework::algorithms
        metada::backends::simple
)

target_link_libraries(particle_filter
    PRIVATE
        metada::traits
        metada::framework::adapters
        metada::base
        metada::framework::runs
        metada::backends::common::utils
        metada::backends::common::observation
        metada::framework::algorithms
        metada::backends::simple
)

target_link_libraries(enkf
    PRIVATE
        metada::traits
        metada::framework::adapters
        metada::base
        metada::framework::runs
        metada::backends::common::utils
        metada::backends::common::observation
        metada::framework::algorithms
        metada::backends::simple
)

target_link_libraries(lwenkf
    PRIVATE
        metada::traits
        metada::framework::adapters
        metada::base
        metada::framework::runs
        metada::backends::common::utils
        metada::backends::common::observation
        metada::framework::algorithms
        metada::backends::simple
)

target_link_libraries(a4denvar
    PRIVATE
        metada::traits
        metada::framework::adapters
        metada::base
        metada::framework::runs
        metada::backends::common::utils
        metada::backends::common::observation
        metada::framework::algorithms
        metada::backends::simple
)

# Set output directory
set_target_properties(etkf
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(letkf
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(particle_filter
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(enkf
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(lwenkf
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(a4denvar
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
# Add format target to ensure consistent code style across the project
AddFormatTarget(etkf ${CMAKE_CURRENT_SOURCE_DIR})
AddFormatTarget(letkf ${CMAKE_CURRENT_SOURCE_DIR})
AddFormatTarget(particle_filter ${CMAKE_CURRENT_SOURCE_DIR})
AddFormatTarget(enkf ${CMAKE_CURRENT_SOURCE_DIR})
AddFormatTarget(lwenkf ${CMAKE_CURRENT_SOURCE_DIR})
AddFormatTarget(a4denvar ${CMAKE_CURRENT_SOURCE_DIR})