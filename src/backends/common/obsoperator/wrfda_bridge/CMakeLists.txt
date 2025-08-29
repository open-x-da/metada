cmake_minimum_required(VERSION 3.20)

project(wrfda_bridge C Fortran)

add_library(wrfda_bridge STATIC)

set_target_properties(wrfda_bridge PROPERTIES
    OUTPUT_NAME wrfda_bridge
)

target_include_directories(wrfda_bridge
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
)

# WRFDA linking configuration
option(WITH_WRFDA "Enable linking real WRFDA obs operator routines" ON)
set(WRFDA_ROOT "" CACHE PATH "Path to WRF/var/da root directory (e.g., D:/linux/WRF/var/da)")

# Allow specifying the prebuilt WRFDA static library directly
set(WRFDA_VAR_LIB "" CACHE FILEPATH "Path to prebuilt libwrfvar.a (e.g., D:/linux/WRF/var/da/build/libwrfvar.a)")

if(WITH_WRFDA)
    # Prefer explicitly provided library path
    if(WRFDA_VAR_LIB)
        if(EXISTS ${WRFDA_VAR_LIB})
            message(STATUS "WRFDA bridge: linking prebuilt library: ${WRFDA_VAR_LIB}")
            target_link_libraries(wrfda_bridge PUBLIC ${WRFDA_VAR_LIB})
        else()
            message(WARNING "WRFDA_VAR_LIB was set but not found: ${WRFDA_VAR_LIB}")
        endif()
    endif()

    # Try to discover library under WRFDA_ROOT if provided
    if(WRFDA_ROOT AND EXISTS ${WRFDA_ROOT})
        find_library(WRFDA_VAR_LIB_DISCOVERED
            NAMES wrfvar libwrfvar.a
            PATHS ${WRFDA_ROOT}/build ${WRFDA_ROOT} ${WRFDA_ROOT}/var/build
            NO_DEFAULT_PATH
        )
        if(WRFDA_VAR_LIB_DISCOVERED)
            message(STATUS "WRFDA bridge: discovered libwrfvar: ${WRFDA_VAR_LIB_DISCOVERED}")
            target_link_libraries(wrfda_bridge PUBLIC ${WRFDA_VAR_LIB_DISCOVERED})
        else()
            message(WARNING "WRFDA bridge: could not find libwrfvar under ${WRFDA_ROOT}. Provide WRFDA_VAR_LIB explicitly.")
        endif()

        # Expose WRFDA includes for .inc files and modules
        target_include_directories(wrfda_bridge PUBLIC ${WRFDA_ROOT} ${WRFDA_ROOT}/var/da)
    endif()

    # Real wrapper sources require WRFDA modules (.mod) to be found
    set(WRFDA_MOD_DIR "" CACHE PATH "Path to WRFDA build directory with .mod files (e.g., D:/linux/WRF/var/da/build)")
    if(WRFDA_MOD_DIR AND EXISTS ${WRFDA_MOD_DIR})
        message(STATUS "WRFDA bridge: using module dir ${WRFDA_MOD_DIR}")
        target_include_directories(wrfda_bridge PUBLIC ${WRFDA_MOD_DIR})
        target_sources(wrfda_bridge PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/wrfda_dispatch.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/wrfda_force_symbols.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/wrfda_bridge.c
        )

        # Some WRFDA builds do not archive da_control.o into libwrfvar.a
        # Link the object explicitly if present to satisfy module globals
        set(WRFDA_DA_CONTROL_OBJ "${WRFDA_MOD_DIR}/da_control.o")
        if (EXISTS ${WRFDA_DA_CONTROL_OBJ})
            message(STATUS "WRFDA bridge: linking extra object ${WRFDA_DA_CONTROL_OBJ}")
            target_link_libraries(wrfda_bridge PUBLIC ${WRFDA_DA_CONTROL_OBJ})
        endif()
    else()
        message(FATAL_ERROR "WRFDA bridge: WRFDA_MOD_DIR must be set to build with real WRFDA observation operators")
    endif()
else()
    # WRFDA is required for this bridge
    message(FATAL_ERROR "WRFDA bridge: WITH_WRFDA must be enabled to build with real WRFDA observation operators")
endif()

add_library(metada::bridges::wrfda ALIAS wrfda_bridge)


