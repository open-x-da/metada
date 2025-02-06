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
  find_program(PERL_EXECUTABLE perl)

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
    # Clean existing coverage data
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_BINARY_DIR}/*.gcda
    COMMAND ${LCOV_COMMAND} --directory ${CMAKE_BINARY_DIR} --zerocounters
    # Run the tests
    COMMAND $<TARGET_FILE:${target}>
    # Capture coverage data
    COMMAND ${LCOV_COMMAND} 
            --directory ${CMAKE_BINARY_DIR}
            --capture 
            --output-file ${target_binary_dir}/coverage.info
            --base-directory ${CMAKE_SOURCE_DIR}
            --no-external
    COMMAND ${LCOV_COMMAND} 
            -r ${target_binary_dir}/coverage.info 
            '/usr/include/*'
            '/mingw64/include/*'
            '${CMAKE_BINARY_DIR}/_deps/*'
            -o ${target_binary_dir}/filtered.info
    COMMAND ${GENHTML_COMMAND} 
            -o ${target_binary_dir}/coverage-${target}
            ${target_binary_dir}/filtered.info 
            --legend
            --base-directory ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove ${target_binary_dir}/coverage.info ${target_binary_dir}/filtered.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating coverage report for ${target}"
  )

  # Add a global coverage target if it doesn't exist
  if(NOT TARGET coverage)
    add_custom_target(coverage)
  endif()
  add_dependencies(coverage coverage-${target})
endfunction()
