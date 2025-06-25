add_library(base 
    DateTime.cpp
    Duration.cpp
)
add_library(metada::base ALIAS base)

target_sources(base
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/NonCopyable.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DateTime.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Duration.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/PointObservation.hpp
)

target_include_directories(base
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Ensure C++20 support
target_compile_features(base PUBLIC cxx_std_20)

AddFormatTarget(base ${CMAKE_CURRENT_SOURCE_DIR})