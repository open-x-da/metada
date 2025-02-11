# State tests
add_executable(state_test
    StateTest.cpp
)

target_link_libraries(state_test
    PRIVATE
        metada::framework::repr
)

metada_add_test_with_coverage(state_test)

AddFormatTarget(state_test ${CMAKE_CURRENT_SOURCE_DIR}) 

# Increment tests
add_executable(increment_test
    IncrementTest.cpp
)

target_link_libraries(increment_test
    PRIVATE
        metada::framework::repr
)

metada_add_test_with_coverage(increment_test)

AddFormatTarget(increment_test ${CMAKE_CURRENT_SOURCE_DIR}) 