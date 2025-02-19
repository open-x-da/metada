include(Coverage)

macro(metada_add_test_with_coverage target)
  AddCoverage(${target})

  target_link_libraries(${target}
    PRIVATE
      GTest::gtest_main
      GTest::gmock)

  gtest_discover_tests(${target} -D TEST_DISCOVERY_TIMEOUT=20)

  #AddMemcheck(${target})
endmacro()