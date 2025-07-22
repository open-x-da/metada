# TL/AD Checks Application
add_executable(tl_ad_checks tl_ad_checks.cpp)

# Variational Application
add_executable(variational variational.cpp)

target_link_libraries(tl_ad_checks
  PRIVATE
      metada::traits
      metada::framework::adapters
      metada::base
      metada::framework::runs
      metada::backends::common::utils
      metada::backends::common::observation
      metada::framework::algorithms
      metada::backends::wrf
)

target_link_libraries(variational
  PRIVATE
      metada::traits
      metada::framework::adapters
      metada::base
      metada::framework::runs
      metada::backends::common::utils
      metada::backends::common::observation
      metada::framework::algorithms
      metada::backends::wrf
)

# Set output directory
set_target_properties(tl_ad_checks
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(variational
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
# Add format target to ensure consistent code style across the project
AddFormatTarget(tl_ad_checks ${CMAKE_CURRENT_SOURCE_DIR}) 
AddFormatTarget(variational ${CMAKE_CURRENT_SOURCE_DIR})