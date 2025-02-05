# Print system information
function(print_system_info)
    print_subheader("Platform Information")
    message(STATUS "  - OS: ${CMAKE_HOST_SYSTEM_NAME} ${CMAKE_HOST_SYSTEM_VERSION}")
    message(STATUS "  - Architecture: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
    
    # Get hostname (platform-specific)
    if(WIN32)
        execute_process(
            COMMAND hostname
            OUTPUT_VARIABLE HOST_NAME
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        execute_process(
            COMMAND cmd /c echo %USERNAME%
            OUTPUT_VARIABLE USER_NAME
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    else()
        execute_process(
            COMMAND uname -n
            OUTPUT_VARIABLE HOST_NAME
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        execute_process(
            COMMAND whoami
            OUTPUT_VARIABLE USER_NAME
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()
    message(STATUS "  - Hostname: ${HOST_NAME}")
    message(STATUS "  - Username: ${USER_NAME}")
endfunction() 