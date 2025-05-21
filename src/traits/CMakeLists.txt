# CMakeLists.txt for traits directory
#
# This file sets up the traits library as a header-only component

add_library(metada_traits INTERFACE)
add_library(metada::traits ALIAS metada_traits)

target_include_directories(metada_traits
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Add compile definitions for the selected backends
if(DEFINED CONFIG_BACKEND)
    string(TOUPPER ${CONFIG_BACKEND} CONFIG_BACKEND_UPPER)
    target_compile_definitions(metada_traits 
        INTERFACE 
        CONFIG_BACKEND_${CONFIG_BACKEND_UPPER}
    )
endif()

if(DEFINED LOGGER_BACKEND)
    string(TOUPPER ${LOGGER_BACKEND} LOGGER_BACKEND_UPPER)
    target_compile_definitions(metada_traits 
        INTERFACE 
        LOGGER_BACKEND_${LOGGER_BACKEND_UPPER}
    )
endif()

# Add any specific compiler flags or definitions needed for traits
target_compile_features(metada_traits
    INTERFACE
        cxx_std_20
) 