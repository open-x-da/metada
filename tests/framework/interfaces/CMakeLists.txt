add_library(interfaces_test INTERFACE)

add_library(metada::framework::interfaces::test
    ALIAS interfaces_test)

# Set up include directories
target_include_directories(interfaces_test
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# List header files for proper IDE indexing and documentation
target_sources(interfaces_test
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/MockModel.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockObservation.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockObsOperator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockTimeStepper.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockAIPredictor.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockBatchProcessor.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockHardwareAccelerator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockGeometry.hpp
)

target_link_libraries(interfaces_test
    INTERFACE
        metada::framework::adapters
        metada::framework::interfaces
        metada::framework::core
)

AddFormatTarget(interfaces_test
    ${CMAKE_CURRENT_SOURCE_DIR})