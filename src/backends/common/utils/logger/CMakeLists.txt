# Logger backend selection
#
# Supports Google Glog or basic console logging.
# Console logger is used as fallback if Glog is not available.

if(Glog_FOUND)
    set(LOGGER_BACKEND "GLOG" CACHE STRING "Select logger backend")
else()
    set(LOGGER_BACKEND "CONSOLE" CACHE STRING "Select logger backend")
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