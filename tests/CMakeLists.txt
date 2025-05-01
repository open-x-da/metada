
if(BUILD_TESTING)
  include(GoogleTest)

  # Add test directories
  add_subdirectory(framework/utils)
  add_subdirectory(framework/adapters)
endif()