include(Coverage)

macro(AddTests target)
  AddCoverage(${target})

  # Set C++17 standard for filesystem support
  # This is required because:
  # 1. std::filesystem is used for managing file and directory operations
  # 2. std::filesystem became standardized in C++17
  # 3. Earlier versions only had experimental::filesystem which required separate linking
  # Setting CXX_STANDARD_REQUIRED ensures the build fails if C++17 is not available
  # rather than silently falling back to an older standard
  set_target_properties(${target} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
  )

  # On older Apple systems (Clang < 9.0), the filesystem library needs to be explicitly linked
  # This is because std::filesystem was experimental before C++17 and required separate linking
  # On newer systems, it's part of the standard library and doesn't need explicit linking
  # Reference: https://en.cppreference.com/w/cpp/filesystem
  if(APPLE AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(${target} PRIVATE c++fs)
  endif()

  target_link_libraries(${target}
    PRIVATE
      GTest::gtest_main
      GTest::gmock)

  gtest_discover_tests(${target})

  #AddMemcheck(${target})
endmacro()