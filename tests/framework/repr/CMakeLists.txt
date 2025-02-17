# State tests
add_executable(state_test
    StateTest.cpp
    EnsembleTest.cpp
)

target_link_libraries(state_test
    PRIVATE
        metada::framework::repr
        metada::backends::repr::mock
)

target_include_directories(state_test
    PRIVATE
        ${CMAKE_SOURCE_DIR}/tests/framework/tools/config
        ${CMAKE_SOURCE_DIR}/tests/framework/tools/logger
)

# Increment tests
#add_executable(increment_test
#       IncrementTest.cpp
#)

#target_link_libraries(increment_test
#    PRIVATE
#        metada::framework::repr
#)

# Model tests
#add_executable(model_test
#    ModelTest.cpp
#)

#target_link_libraries(model_test
#    PRIVATE
#        metada::framework::repr
#)

metada_add_test_with_coverage(state_test)
#metada_add_test_with_coverage(increment_test)
#metada_add_test_with_coverage(model_test)

AddFormatTarget(state_test ${CMAKE_CURRENT_SOURCE_DIR}) 
#AddFormatTarget(increment_test ${CMAKE_CURRENT_SOURCE_DIR}) 
#AddFormatTarget(model_test ${CMAKE_CURRENT_SOURCE_DIR}) 