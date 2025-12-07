include(TestingAndCoverage)

# Define a helper function to set up test executables with C++20 support
function(setup_test_executable target source)
    add_executable(${target} ${source})
    
    target_link_libraries(${target}
        PRIVATE
            metada::traits
            metada::framework::adapters
            metada::framework::runs
            metada::base
            metada::framework::algorithms
            metada::backends::gmock
    ) 
  
    # Add test with coverage
    metada_add_test_with_coverage(${target})
endfunction()

# Define a function to set up MPI test execution
# This should be called after setup_test_executable for MPI tests
# Usage: setup_mpi_test(target [num_procs])
function(setup_mpi_test target)
    # Parse optional num_procs argument (default to 4)
    set(num_procs 4)
    if(ARGC GREATER 1)
        set(num_procs ${ARGV1})
    endif()
    
    # This function is now a wrapper that calls the MPI-aware test setup
    # The actual MPI wrapping happens in wrap_mpi_tests() which should be
    # called at the end of the CMakeLists.txt
    message(STATUS "MPI test setup for ${target} (${num_procs} processes)")
endfunction()

# Define a function to add formatting for test files
function(add_test_format_target target directory)
    AddFormatTarget(${target} ${directory})
endfunction() 