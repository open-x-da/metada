include(Coverage)

macro(AddTests target)
  AddCoverage(${target})

  target_link_libraries(${target}
    PRIVATE
      GTest::gtest_main
      GTest::gmock)

  gtest_discover_tests(${target})

  #AddMemcheck(${target})
endmacro()