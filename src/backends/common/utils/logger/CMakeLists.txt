# Logger backend selection
#
# Supports ng-log or basic console logging.
# Console logger is used as fallback if ng-log is not available.

if(NOT DEFINED LOGGER_BACKEND)
    if(ng-log_FOUND)
        set(LOGGER_BACKEND "NGLOG" CACHE STRING "Select logger backend (NGLOG or CONSOLE)")
        set_property(CACHE LOGGER_BACKEND PROPERTY STRINGS "NGLOG" "CONSOLE")
    else()
        set(LOGGER_BACKEND "CONSOLE" CACHE STRING "Select logger backend")
    endif()
else()
    # Validate user-selected backend
    string(TOUPPER ${LOGGER_BACKEND} LOGGER_BACKEND_UPPER)
    if(LOGGER_BACKEND_UPPER STREQUAL "NGLOG" AND NOT ng-log_FOUND)
        message(FATAL_ERROR "NGLOG logger backend requested but nglog is not available")
    elseif(NOT LOGGER_BACKEND_UPPER STREQUAL "NGLOG" AND NOT LOGGER_BACKEND_UPPER STREQUAL "CONSOLE")
        message(FATAL_ERROR "Invalid logger backend: ${LOGGER_BACKEND}. Valid options are: NGLOG, CONSOLE")
    endif()
    set(LOGGER_BACKEND ${LOGGER_BACKEND_UPPER})
endif()

# Create an internal logger interface target
add_library(logger_backend INTERFACE)

# Add the selected backend as a dependency
string(TOLOWER ${LOGGER_BACKEND} LOGGER_BACKEND_LOWER)
add_subdirectory(${LOGGER_BACKEND_LOWER})

# Link the selected backend to the interface library
target_link_libraries(logger_backend
    INTERFACE
        logger_${LOGGER_BACKEND_LOWER}
)