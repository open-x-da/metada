add_subdirectory(common)
add_subdirectory(lorenz63)
add_subdirectory(gmock)
add_subdirectory(lite)
add_subdirectory(simple)

# WRF backend requires WRFDA support
if(WITH_WRFDA)
    add_subdirectory(wrf)
endif()