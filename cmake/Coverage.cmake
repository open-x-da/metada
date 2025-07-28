# This file is adopted from Modern CMake for C++, published by Packt.

# Enables code coverage instrumentation for a target in Debug mode
# Adds compiler-specific coverage flags and disables inlining for accurate coverage data
function(enable_coverage target)
  # Skip if lcov is not found
  if(NOT Lcov_FOUND)
    return()
  endif()

  if(NOT CMAKE_BUILD_TYPE STREQUAL Debug)
    return()
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # GCC and Clang use --coverage flag
    target_compile_options(${target} PRIVATE 
      --coverage
      -fno-inline
    )
    target_link_options(${target} PUBLIC --coverage)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    # Updated Intel oneAPI coverage flags
    target_compile_options(${target} PRIVATE
      -fprofile-instr-generate
      -fcoverage-mapping
      -fno-inline
    )
    target_link_options(${target} PUBLIC
      -fprofile-instr-generate
      -fcoverage-mapping
    )
  else()
    message(WARNING "Code coverage not supported for compiler: ${CMAKE_CXX_COMPILER_ID}")
  endif()
endfunction()

# Creates a coverage report target that:
# 1. Resets coverage counters
# 2. Runs the tests
# 3. Captures coverage data
# 4. Filters out system headers
# 5. Generates HTML report
function(AddCoverage target)
  # Check if lcov and genhtml are available
  if(NOT Lcov_FOUND)
    message(STATUS "Skipping coverage for ${target} - lcov/genhtml not found")
    return()
  endif()
  
  if(WIN32 AND PERL_EXECUTABLE)
    set(LCOV_COMMAND ${CMAKE_COMMAND} -E env "LC_ALL=C" "${PERL_EXECUTABLE}" "${LCOV_PATH}")
    set(GENHTML_COMMAND ${CMAKE_COMMAND} -E env "LC_ALL=C" "${PERL_EXECUTABLE}" "${GENHTML_PATH}")
  else()
    set(LCOV_COMMAND "${LCOV_PATH}")
    set(GENHTML_COMMAND "${GENHTML_PATH}")
  endif()

  # Enable coverage for the target
  enable_coverage(${target})

  # Get the binary directory for this target
  get_target_property(target_binary_dir ${target} BINARY_DIR)
  if(NOT target_binary_dir)
    set(target_binary_dir ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  add_custom_target(
    coverage-${target}
    # Clean up the old gcda
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} 
            -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR} 
            -P ${CMAKE_SOURCE_DIR}/cmake/scripts/clean_gcda.cmake
    # Capture coverage data
    COMMAND ${LCOV_COMMAND} 
            --directory .
            --capture 
            --output-file coverage.info
            --base-directory ${CMAKE_SOURCE_DIR}
    # Filter out system headers
    COMMAND ${LCOV_COMMAND} 
            -r coverage.info 
            "C:/msys64/mingw64/*"
            "*/mingw64/*"
            "*/include/c++/*"
            "*/bits/*"
            "*/gtest/*"
            "*/gmock/*"
            "*/usr/include/*" 
            "*/usr/lib/*" 
            "*/opt/*" 
            "${CMAKE_BINARY_DIR}/_deps/*"
            -o filtered.info
    # Generate HTML report
    COMMAND ${GENHTML_COMMAND} 
            --output-directory ${CMAKE_BINARY_DIR}/coverage/${target}
            --prefix ${CMAKE_SOURCE_DIR}
            --rc geninfo_unexecuted_blocks=1 
            --ignore-errors mismatch
            filtered.info
            --legend
    # Clean up
    COMMAND ${CMAKE_COMMAND} -E remove coverage.info filtered.info
    WORKING_DIRECTORY ${target_binary_dir}
    COMMENT "Generating coverage report for ${target}"
  )

  # Add a global coverage target if it doesn't exist
  if(NOT TARGET coverage)
    add_custom_target(coverage)
  endif()
  add_dependencies(coverage coverage-${target})
endfunction()
