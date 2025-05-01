add_library(utils 
    DateTime.cpp
    Duration.cpp
)
add_library(metada::framework::utils ALIAS utils)

target_sources(utils
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/NonCopyable.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DateTime.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Duration.hpp
)

target_include_directories(utils
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Ensure C++20 support
target_compile_features(utils PUBLIC cxx_std_20)

AddFormatTarget(utils ${CMAKE_CURRENT_SOURCE_DIR})