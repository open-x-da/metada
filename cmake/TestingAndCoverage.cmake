include(Coverage)

macro(metada_add_test_with_coverage target)
  AddCoverage(${target})

  target_link_libraries(${target}
    PRIVATE
      GTest::gtest_main
      GTest::gmock)

  # Add a pre-test hook to clean .gcda files before running this test
  add_custom_command(
    TARGET ${target}
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR} -P ${CMAKE_SOURCE_DIR}/cmake/scripts/clean_gcda.cmake
    COMMENT "Cleaning coverage data before building ${target}"
  )

  gtest_discover_tests(${target} -D TEST_DISCOVERY_TIMEOUT=60)

  #AddMemcheck(${target})
endmacro()