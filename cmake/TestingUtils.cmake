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
            metada::backends::gmock
    ) 
  
    # Add test with coverage
    metada_add_test_with_coverage(${target})
endfunction()

# Define a function to add formatting for test files
function(add_test_format_target target directory)
    AddFormatTarget(${target} ${directory})
endfunction() 