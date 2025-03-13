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
        ${CMAKE_CURRENT_SOURCE_DIR}/MockState.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockModel.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockObservation.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockObsOperator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockTimeStepper.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockAIPredictor.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockBatchProcessor.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockHardwareAccelerator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MockGeometry.hpp
        # C++20 Module mock implementations
        ${CMAKE_CURRENT_SOURCE_DIR}/MockGeometry.ixx
)

# Set module file properties
set_module_file_properties("${CMAKE_CURRENT_SOURCE_DIR}/MockGeometry.ixx")

# Set C++20 properties for the target
target_compile_features(interfaces_test INTERFACE cxx_std_20)

# For INTERFACE libraries, we need different property setting approach
set_property(TARGET interfaces_test PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_20)

# Add compiler-specific flags for C++20 modules
add_module_compiler_flags(interfaces_test)

# Add compiler-specific flags for C++20 concepts via generator expressions
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(interfaces_test INTERFACE 
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<VERSION_LESS:${CMAKE_CXX_COMPILER_VERSION},10.0>>:-fconcepts>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(interfaces_test INTERFACE 
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<VERSION_LESS:${CMAKE_CXX_COMPILER_VERSION},10.0>>:-fconcepts-ts>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(interfaces_test INTERFACE /Zc:__cplusplus /permissive-)
endif()

target_link_libraries(interfaces_test
    INTERFACE
        metada::framework::common
)

AddFormatTarget(interfaces_test
    ${CMAKE_CURRENT_SOURCE_DIR})