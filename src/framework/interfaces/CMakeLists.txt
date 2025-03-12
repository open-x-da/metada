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
        ${CMAKE_CURRENT_SOURCE_DIR}/GeometryIterator.hpp
)

# Set C++20 properties for the target
target_compile_features(interfaces INTERFACE cxx_std_20)

# For INTERFACE libraries, we need different property setting approach
set_property(TARGET interfaces PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_20)

# Add compiler-specific flags for C++20 concepts via generator expressions
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(interfaces INTERFACE 
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<VERSION_LESS:${CMAKE_CXX_COMPILER_VERSION},10.0>>:-fconcepts>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(interfaces INTERFACE 
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<VERSION_LESS:${CMAKE_CXX_COMPILER_VERSION},10.0>>:-fconcepts-ts>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(interfaces INTERFACE /Zc:__cplusplus /permissive-)
endif()

AddFormatTarget(interfaces ${CMAKE_CURRENT_SOURCE_DIR})