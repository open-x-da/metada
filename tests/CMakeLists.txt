
if(BUILD_TESTING)
  include(GoogleTest)

  # Add test directories
  add_subdirectory(framework/core)
  add_subdirectory(framework/adapters)
endif()