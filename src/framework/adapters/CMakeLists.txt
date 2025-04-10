# Define the scientific representation library as interface (header-only)
add_library(adapters INTERFACE)

# Create alias target
add_library(metada::framework::adapters ALIAS adapters)

# Set up include directories for interface library
target_include_directories(adapters
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common/utils/config>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common/utils/logger>
        )

# List header files for proper IDE indexing and documentation
target_sources(adapters
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/State.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/StateConcepts.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Increment.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Observation.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ObsOperator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Ensemble.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Model.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Geometry.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/GeometryPointIterator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/GeometryIterator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/utils/config/Config.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/utils/config/ConfigConcepts.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/utils/config/ConfigValue.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/utils/logger/Logger.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/utils/logger/LogLevel.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/utils/logger/LoggerConcepts.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/utils/logger/LogStream.hpp
)

AddFormatTarget(adapters ${CMAKE_CURRENT_SOURCE_DIR})