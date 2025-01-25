# Print system information
function(print_system_info)
    message(STATUS "System Information:")
    message(STATUS "  - OS: ${CMAKE_HOST_SYSTEM_NAME} ${CMAKE_HOST_SYSTEM_VERSION}")
    message(STATUS "  - Architecture: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
    message(STATUS "  - CMake Version: ${CMAKE_VERSION}")
    message(STATUS "  - CMake Path: ${CMAKE_COMMAND}")
endfunction() 