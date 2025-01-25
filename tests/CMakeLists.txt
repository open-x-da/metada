# Find required packages
metada_find_package(GTest)

# Add test categories
add_subdirectory(cuda)
add_subdirectory(logging)
add_subdirectory(python)