# Define the scientific representation library as interface (header-only)
add_library(adapters INTERFACE)

# Create alias target
add_library(metada::framework::adapters ALIAS adapters)

# Set up include directories for interface library
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
        ${CMAKE_CURRENT_SOURCE_DIR}/Geometry.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/GeometryPointIterator.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/model/TimeStepperImpl.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/model/AIPredictorImpl.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/model/BatchProcessorImpl.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/model/HardwareAcceleratorImpl.hpp
)

# Set C++20 properties for the target
target_compile_features(adapters INTERFACE cxx_std_20)

# Add compiler-specific flags for C++20 concepts via generator expressions
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(adapters INTERFACE 
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<VERSION_LESS:${CMAKE_CXX_COMPILER_VERSION},10.0>>:-fconcepts>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(adapters INTERFACE 
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<VERSION_LESS:${CMAKE_CXX_COMPILER_VERSION},10.0>>:-fconcepts-ts>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(adapters INTERFACE /Zc:__cplusplus /permissive-)
endif()

# Ensure we depend on interfaces which also has C++20 configured
target_link_libraries(adapters INTERFACE metada::framework::interfaces)

AddFormatTarget(adapters ${CMAKE_CURRENT_SOURCE_DIR})