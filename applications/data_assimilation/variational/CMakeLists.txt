# TL/AD Checks Application
add_executable(tl_ad_checks tl_ad_checks.cpp)

# Variational Application
add_executable(variational variational.cpp)

# Incremental Variational Application
add_executable(incremental_variational incremental_variational.cpp)

target_link_libraries(tl_ad_checks
  PRIVATE
      metada::traits
      metada::framework::adapters
      metada::base
      metada::framework::runs
      metada::backends::common::utils
      metada::backends::common::io
      metada::backends::common::observation
      metada::backends::common::obsoperator
      metada::bridges::wrfda
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
      metada::backends::common::io
      metada::backends::common::obsoperator
      metada::bridges::wrfda
      metada::framework::algorithms
      metada::backends::wrf
)

target_link_libraries(incremental_variational
  PRIVATE
      metada::traits
      metada::framework::adapters
      metada::base
      metada::framework::runs
      metada::backends::common::utils
      metada::backends::common::observation
      metada::backends::common::io
      metada::backends::common::obsoperator
      metada::bridges::wrfda
      metada::framework::algorithms
      metada::backends::wrf
)

# If WRFDA is enabled, use linker group to resolve circular deps between wrfda_bridge and libwrfvar
if(WITH_WRFDA AND WRFDA_VAR_LIB)
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Explicitly add da_control.o if present and not archived in libwrfvar
    set(WRFDA_DA_CONTROL_OBJ "${WRFDA_MOD_DIR}/da_control.o")
    if (EXISTS ${WRFDA_DA_CONTROL_OBJ})
      target_link_libraries(variational PRIVATE ${WRFDA_DA_CONTROL_OBJ})
      target_link_libraries(incremental_variational PRIVATE ${WRFDA_DA_CONTROL_OBJ})
    endif()
    target_link_options(variational PRIVATE -Wl,--start-group)
    target_link_libraries(variational PRIVATE ${WRFDA_VAR_LIB})
    target_link_options(variational PRIVATE -Wl,--end-group)
    target_link_options(incremental_variational PRIVATE -Wl,--start-group)
    target_link_libraries(incremental_variational PRIVATE ${WRFDA_VAR_LIB})
    target_link_options(incremental_variational PRIVATE -Wl,--end-group)
  endif()
endif()

# Set output directory
set_target_properties(tl_ad_checks
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(variational
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(incremental_variational
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
# Add format target to ensure consistent code style across the project
AddFormatTarget(tl_ad_checks ${CMAKE_CURRENT_SOURCE_DIR}) 
AddFormatTarget(variational ${CMAKE_CURRENT_SOURCE_DIR})
AddFormatTarget(incremental_variational ${CMAKE_CURRENT_SOURCE_DIR})