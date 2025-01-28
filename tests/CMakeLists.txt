# Add test categories
if(WIN32 AND MINGW)
    add_subdirectory(cuda)
    add_subdirectory(logging)
    # Skip Python tests on Windows with MSYS2
    message(STATUS "Skipping Python tests on Windows with MSYS2")
else()
    add_subdirectory(cuda)
    add_subdirectory(logging)
    add_subdirectory(python)
endif()