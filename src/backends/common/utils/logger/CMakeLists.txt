# Define available logger backends and set default
if(glog_FOUND)
    set(LOGGER_BACKEND "GLOG" CACHE STRING "Select logger backend")
else()
    set(LOGGER_BACKEND "CONSOLE" CACHE STRING "Select logger backend") 
endif()

# Set available options based on found packages
set(LOGGER_BACKENDS "")
if(glog_FOUND)
    list(APPEND LOGGER_BACKENDS "GLOG")
endif()
list(APPEND LOGGER_BACKENDS "CONSOLE")

# Set the allowed values for LOGGER_BACKEND
set_property(CACHE LOGGER_BACKEND PROPERTY STRINGS ${LOGGER_BACKENDS})

# Add selected backend subdirectory
if(LOGGER_BACKEND IN_LIST LOGGER_BACKENDS)
    string(TOLOWER ${LOGGER_BACKEND} LOGGER_BACKEND_LOWER)
    add_subdirectory(${LOGGER_BACKEND_LOWER})
else()
    message(FATAL_ERROR "Invalid LOGGER_BACKEND value: ${LOGGER_BACKEND}. Must be one of: ${LOGGER_BACKENDS}")
endif()

# Create backend logger library
add_library(logger_backends INTERFACE)
add_library(metada::backends::common::utils::logger ALIAS logger_backends)

target_include_directories(logger_backends
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Configure compile definitions and link libraries based on selected backend
target_compile_definitions(logger_backends
    INTERFACE
        USE_${LOGGER_BACKEND}_LOGGER
)

target_link_libraries(logger_backends
    INTERFACE
        metada::backends::common::utils::logger::${LOGGER_BACKEND_LOWER}
)