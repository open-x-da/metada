# LETKF tests
setup_test_executable(etkf_test ETKFTest.cpp)
add_test_format_target(etkf_test ${CMAKE_CURRENT_SOURCE_DIR}) 

# Metrics tests
setup_test_executable(metrics_test MetricsTest.cpp)
add_test_format_target(metrics_test ${CMAKE_CURRENT_SOURCE_DIR})