include(Coverage)

macro(metada_add_test_with_coverage target)
  # Add coverage if lcov is available
  if(Lcov_FOUND)
    AddCoverage(${target})
    
    # Add a pre-test hook to clean .gcda files before running this test
    add_custom_command(
      TARGET ${target}
      PRE_BUILD
      COMMAND ${CMAKE_COMMAND} 
              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} 
              -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR} 
              -P ${CMAKE_SOURCE_DIR}/cmake/scripts/clean_gcda.cmake
      COMMENT "Cleaning coverage data before building ${target}"
    )
  else()
    message(STATUS "Coverage analysis disabled for ${target} - lcov/genhtml not found")
  endif()

  target_link_libraries(${target}
    PRIVATE
      GTest::gtest_main
      GTest::gmock
  )

  gtest_discover_tests(
    ${target}
    PROPERTIES
      TIMEOUT 300
    DISCOVERY_TIMEOUT 300
  )

  # AddMemcheck(${target})
endmacro()

# Macro to add test with coverage and MPI support
# This wraps the test command with mpirun for parallel execution
# Note: The executable must already be created before calling this macro
macro(metada_add_test_with_coverage_mpi target)
  # Parse optional num_procs argument (default to 4)
  set(num_procs 4)
  if(ARGC GREATER 1)
    set(num_procs ${ARGV1})
  endif()

  # Ensure target exists
  if(NOT TARGET ${target})
    message(FATAL_ERROR "Target ${target} must be created before calling metada_add_test_with_coverage_mpi")
  endif()

  # Add coverage if lcov is available
  if(Lcov_FOUND)
    AddCoverage(${target})
    
    # Add a pre-test hook to clean .gcda files before running this test
    add_custom_command(
      TARGET ${target}
      PRE_BUILD
      COMMAND ${CMAKE_COMMAND} 
              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} 
              -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR} 
              -P ${CMAKE_SOURCE_DIR}/cmake/scripts/clean_gcda.cmake
      COMMENT "Cleaning coverage data before building ${target}"
    )
  else()
    message(STATUS "Coverage analysis disabled for ${target} - lcov/genhtml not found")
  endif()

  target_link_libraries(${target}
    PRIVATE
      GTest::gtest_main
      GTest::gmock
  )

  # Find MPI launcher if not already set
  if(NOT MPIEXEC)
    find_program(MPIEXEC
      NAMES mpirun mpiexec
      PATHS ${MPI_HOME}/bin
            C:/msys64/mingw64/bin
            C:/msys64/clang64/bin
            /usr/bin
            /usr/local/bin
      DOC "MPI launcher for running parallel tests"
    )
  endif()

  if(MPIEXEC)
    # Set number of processes flag if not set
    if(NOT MPIEXEC_NUMPROC_FLAGS)
      set(MPIEXEC_NUMPROC_FLAGS "-n")
    endif()

    # Use gtest_discover_tests but wrap the executable with mpirun
    # by setting CROSSCOMPILING_EMULATOR on the target
    # This should work, but if it doesn't, we'll use a post-build script
    set_target_properties(${target} PROPERTIES
      CROSSCOMPILING_EMULATOR "${MPIEXEC};${MPIEXEC_NUMPROC_FLAGS};${num_procs}"
    )

    # Discover tests - CROSSCOMPILING_EMULATOR should prepend mpirun
    gtest_discover_tests(
      ${target}
      TEST_PREFIX "${target}."
      PROPERTIES
        TIMEOUT 300
      DISCOVERY_TIMEOUT 300
    )

    # Verify that CROSSCOMPILING_EMULATOR is set
    get_target_property(emulator ${target} CROSSCOMPILING_EMULATOR)
    if(emulator)
      message(STATUS "  CROSSCOMPILING_EMULATOR set: ${emulator}")
    else()
      message(WARNING "  CROSSCOMPILING_EMULATOR not set on ${target}")
    endif()

    # Also try to wrap tests directly as a fallback
    # Get all tests and wrap them
    get_property(test_list DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
    
    foreach(test_name IN LISTS test_list)
      if(test_name MATCHES "^${target}\\.")
        get_test_property(${test_name} COMMAND test_cmd)
        if(test_cmd)
          # Check if already wrapped
          list(GET test_cmd 0 first_arg)
          if(NOT first_arg STREQUAL ${MPIEXEC})
            # Not wrapped yet, wrap it
            set(new_cmd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAGS} ${num_procs} ${test_cmd})
            set_tests_properties(${test_name} PROPERTIES COMMAND "${new_cmd}")
            message(STATUS "  Manually wrapped: ${test_name}")
          endif()
        endif()
      endif()
    endforeach()
    
    message(STATUS "MPI test configuration for ${target}:")
    message(STATUS "  MPI launcher: ${MPIEXEC}")
    message(STATUS "  Number of processes: ${num_procs}")
    message(STATUS "  Process flag: ${MPIEXEC_NUMPROC_FLAGS}")
    message(STATUS "  Tests will be wrapped with: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAGS} ${num_procs}")
  else()
    message(WARNING "MPI launcher not found. ${target} will run serially.")
    # Fall back to regular test discovery
    metada_add_test_with_coverage(${target})
  endif()
endmacro()

# Function to wrap MPI test commands after discovery
# This should be called at the end of the CMakeLists.txt that defines MPI tests
function(wrap_mpi_tests target)
  if(NOT MPI_TEST_CONFIG_${target}_LAUNCHER)
    return()
  endif()

  set(mpi_launcher ${MPI_TEST_CONFIG_${target}_LAUNCHER})
  set(mpi_flags ${MPI_TEST_CONFIG_${target}_FLAGS})
  set(num_procs ${MPI_TEST_CONFIG_${target}_NUMPROC})

  # Get all tests and wrap those that match the target prefix
  get_property(test_list DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
  
  foreach(test_name IN LISTS test_list)
    if(test_name MATCHES "^${target}\\.")
      # Get the current test command
      get_test_property(${test_name} COMMAND test_cmd)
      if(test_cmd)
        # Wrap the command with mpirun
        # test_cmd is a list, so we need to preserve it properly
        set(new_cmd ${mpi_launcher} ${mpi_flags} ${num_procs} ${test_cmd})
        set_tests_properties(${test_name} PROPERTIES COMMAND "${new_cmd}")
        message(STATUS "Wrapped MPI test: ${test_name}")
      endif()
    endif()
  endforeach()
endfunction()