# CMakeLists.txt for traits directory
#
# This file sets up the traits library as a header-only component

add_library(metada_traits INTERFACE)
add_library(metada::traits ALIAS metada_traits)

target_include_directories(metada_traits
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Add any specific compiler flags or definitions needed for traits
target_compile_features(metada_traits
    INTERFACE
        cxx_std_20
) 