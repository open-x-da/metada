metada_find_package(bufr CONFIG)

# Create the common io library (only public interface)
add_library(common_io INTERFACE)
add_library(metada::backends::common::io ALIAS common_io)

# Add include directories for utility headers
target_include_directories(common_io
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Link with required libraries
target_link_libraries(common_io
    INTERFACE
        bufr::bufr_4
)

# Add format target
AddFormatTarget(common_io ${CMAKE_CURRENT_SOURCE_DIR})
