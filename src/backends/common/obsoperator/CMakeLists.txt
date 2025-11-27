# Create common obsoperator library
add_library(obsoperator INTERFACE)
add_library(metada::backends::common::obsoperator ALIAS obsoperator)

# Set include directories
target_include_directories(obsoperator
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Link dependencies
target_link_libraries(obsoperator
    INTERFACE
        metada::framework::adapters)

# Add format target
AddFormatTarget(obsoperator ${CMAKE_CURRENT_SOURCE_DIR}) 

# Build the WRFDA bridge (requires WRFDA support)
if(WITH_WRFDA)
    add_subdirectory(wrfda_bridge)
endif()