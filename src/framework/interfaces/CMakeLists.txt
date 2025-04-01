# Define the scientific representation library
add_library(interfaces INTERFACE)

# Create alias target
add_library(metada::framework::interfaces ALIAS interfaces)

# Set up include directories
target_include_directories(interfaces
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# List header files for proper IDE indexing and documentation
target_sources(interfaces
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/IState.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/IObservation.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/IObsOperator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/IModel.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ITimeStepper.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/IAIPredictor.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/IBatchProcessor.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/IHardwareAccelerator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/IGeometry.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/IGeometryIterator.hpp
)

AddFormatTarget(interfaces ${CMAKE_CURRENT_SOURCE_DIR})