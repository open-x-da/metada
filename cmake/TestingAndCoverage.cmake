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