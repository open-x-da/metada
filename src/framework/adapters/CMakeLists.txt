# Define the scientific representation library
add_library(adapters INTERFACE)

# Create alias target
add_library(metada::framework::adapters ALIAS adapters)

# Set up include directories
target_include_directories(adapters
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# List header files for proper IDE indexing and documentation
target_sources(adapters
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/State.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Increment.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Observation.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ObsOperator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Ensemble.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Model.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/model/TimeStepperImpl.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/model/AIPredictorImpl.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/model/BatchProcessorImpl.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/model/HardwareAcceleratorImpl.hpp
)

AddFormatTarget(adapters ${CMAKE_CURRENT_SOURCE_DIR})