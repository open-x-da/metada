add_library(base 
    DateTime.cpp
    Duration.cpp
)
add_library(metada::base ALIAS base)

target_include_directories(base
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Ensure C++20 support
target_compile_features(base PUBLIC cxx_std_20)

AddFormatTarget(base ${CMAKE_CURRENT_SOURCE_DIR})