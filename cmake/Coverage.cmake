# This file is adopted from Modern CMake for C++, published by Packt.

# Enables code coverage instrumentation for a target in Debug mode
# Adds coverage flags and disables inlining for accurate coverage data
function(enable_coverage target)
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    target_compile_options(${target} PRIVATE 
      --coverage
      -fno-inline
    )
    target_link_options(${target} PUBLIC --coverage)
  endif()
endfunction()

# Cleans up coverage data files (*.gcda) before building the target
# This ensures clean coverage data for each test run
function(CleanCoverage target)
  add_custom_target(clean_coverage_${target}
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning coverage data for ${target} in ${CMAKE_CURRENT_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/scripts/clean_coverage.cmake
  )
  
  add_dependencies(${target} clean_coverage_${target})
endfunction()

# Creates a coverage report target that:
# 1. Resets coverage counters
# 2. Runs the tests
# 3. Captures coverage data
# 4. Filters out system headers
# 5. Generates HTML report
function(AddCoverage target)
  find_program(LCOV_PATH lcov REQUIRED)
  find_program(GENHTML_PATH genhtml REQUIRED)
  add_custom_target(
    coverage-${target}
    COMMAND ${LCOV_PATH} -d . --zerocounters
    COMMAND $<TARGET_FILE:${target}>
    COMMAND ${LCOV_PATH} -d . --capture -o coverage.info
    COMMAND ${LCOV_PATH} 
            -r coverage.info 
            '/usr/include/*'
            -o filtered.info
    COMMAND ${GENHTML_PATH} 
            -o coverage-${target}
            filtered.info 
            --legend
    COMMAND rm -rf coverage.info filtered.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endfunction()
