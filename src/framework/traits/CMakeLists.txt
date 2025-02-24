add_library(traits INTERFACE)
add_library(metada::framework::traits ALIAS traits)

target_include_directories(traits
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)