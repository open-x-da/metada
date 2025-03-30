add_library(core INTERFACE)
add_library(metada::framework::core ALIAS core)

target_include_directories(core
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)