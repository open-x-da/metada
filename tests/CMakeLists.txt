
if(BUILD_TESTING)
  include(GoogleTest)

  # Add test directories
  add_subdirectory(framework/base)
  add_subdirectory(framework/adapters)
endif()