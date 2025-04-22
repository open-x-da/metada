add_library(core 
    DateTime.cpp
)
add_library(metada::framework::core ALIAS core)

target_sources(core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/NonCopyable.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DateTime.hpp
)

target_include_directories(core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Ensure C++20 support
target_compile_features(core PUBLIC cxx_std_20)

AddFormatTarget(core ${CMAKE_CURRENT_SOURCE_DIR})