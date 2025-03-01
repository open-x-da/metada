# State tests
add_executable(state_test
    StateTest.cpp
)

target_link_libraries(state_test
    PRIVATE
        metada::framework::interfaces
        metada::framework::adapters
        metada::framework::interfaces::test
        metada::framework::runs
        metada::framework::traits
)

target_include_directories(state_test
    PRIVATE
        ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/config
        ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/logger
)

# Ensemble tests
add_executable(ensemble_test
    EnsembleTest.cpp
)

target_link_libraries(ensemble_test
    PRIVATE
        metada::framework::interfaces
        metada::framework::adapters
        metada::framework::interfaces::test
)

target_include_directories(ensemble_test
    PRIVATE
        ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/config
        ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/logger
)

# Increment tests
add_executable(increment_test
    IncrementTest.cpp
)

target_link_libraries(increment_test
    PRIVATE
        metada::framework::interfaces
        metada::framework::adapters
        metada::framework::interfaces::test
)

target_include_directories(increment_test
    PRIVATE
        ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/config
        ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/logger
)

# Model tests
#add_executable(model_test
#    ModelTest.cpp
#)

#target_link_libraries(model_test
#    PRIVATE
#        metada::framework::repr
#)

# IObservations tests
add_executable(observation_test
    ObservationTest.cpp
)

target_link_libraries(observation_test
    PRIVATE
        metada::framework::interfaces
        metada::framework::adapters
        metada::framework::interfaces::test
)

target_include_directories(observation_test
    PRIVATE
        ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/config
        ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/logger
)

# IObsOperator tests
#add_executable(obsoperator_test
#    ObsOperatorTest.cpp
#)

#target_link_libraries(obsoperator_test
#    PRIVATE
#        metada::framework::interfaces
#        metada::framework::adapters
#        metada::framework::interfaces::test
#)

#target_include_directories(obsoperator_test
#    PRIVATE
#        ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/config
#        ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/logger
#)

metada_add_test_with_coverage(state_test)
metada_add_test_with_coverage(ensemble_test)
metada_add_test_with_coverage(increment_test)
#metada_add_test_with_coverage(model_test)
metada_add_test_with_coverage(observation_test)
#metada_add_test_with_coverage(obsoperator_test)

# AddFormatTarget(state_test ${CMAKE_CURRENT_SOURCE_DIR})
# AddFormatTarget(ensemble_test ${CMAKE_CURRENT_SOURCE_DIR})
# AddFormatTarget(increment_test ${CMAKE_CURRENT_SOURCE_DIR})
# AddFormatTarget(model_test ${CMAKE_CURRENT_SOURCE_DIR})
AddFormatTarget(observation_test ${CMAKE_CURRENT_SOURCE_DIR})
# AddFormatTarget(obsoperator_test ${CMAKE_CURRENT_SOURCE_DIR})