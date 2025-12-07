add_library(base 
    DateTime.cpp
    Duration.cpp
    MPIContext.cpp
)
add_library(metada::base ALIAS base)

target_include_directories(base
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Ensure C++20 support
target_compile_features(base PUBLIC cxx_std_20)

# Link MPI if available
if(METADA_USE_MPI)
    target_link_libraries(base PUBLIC MPI::MPI_CXX)
endif()

AddFormatTarget(base ${CMAKE_CURRENT_SOURCE_DIR})