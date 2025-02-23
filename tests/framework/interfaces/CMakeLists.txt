add_library(interfaces_test INTERFACE)

add_library(metada::framework::interfaces::test
    ALIAS interfaces_test)

# Set up include directories
target_include_directories(interfaces_test
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(interfaces_test
    INTERFACE
        metada::framework::common
)

AddFormatTarget(interfaces_test
    ${CMAKE_CURRENT_SOURCE_DIR})