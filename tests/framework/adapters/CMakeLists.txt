# Define a helper function to set up test executables with C++20 support
function(setup_test_executable target source)
    add_executable(${target} ${source})
    
    target_link_libraries(${target}
        PRIVATE
            metada::traits
            metada::framework::interfaces
            metada::framework::adapters
            metada::framework::interfaces::test
            metada::framework::runs
            metada::framework::traits
            metada::framework::common::utils::logger
    )
    
    target_include_directories(${target}
        PRIVATE
            ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/config
            ${CMAKE_SOURCE_DIR}/tests/framework/common/utils/logger
    )
    
  
    # Add test with coverage
    metada_add_test_with_coverage(${target})
endfunction()

# Config tests
setup_test_executable(config_test ConfigTest.cpp)

# Logger tests
setup_test_executable(logger_test LoggerTest.cpp)

# State tests
#setup_test_executable(state_test StateTest.cpp)

# Ensemble tests
#setup_test_executable(ensemble_test EnsembleTest.cpp)

# Increment tests
#setup_test_executable(increment_test IncrementTest.cpp)

# Model tests
#setup_test_executable(model_test ModelTest.cpp)

# IObservations tests
#setup_test_executable(observation_test ObservationTest.cpp)

# IObsOperator tests
#setup_test_executable(obsoperator_test ObsOperatorTest.cpp)

# Geometry tests
#setup_test_executable(geometry_test GeometryTest.cpp)

AddFormatTarget(config_test ${CMAKE_CURRENT_SOURCE_DIR})
AddFormatTarget(logger_test ${CMAKE_CURRENT_SOURCE_DIR})
# AddFormatTarget(state_test ${CMAKE_CURRENT_SOURCE_DIR})
# AddFormatTarget(ensemble_test ${CMAKE_CURRENT_SOURCE_DIR})
# AddFormatTarget(increment_test ${CMAKE_CURRENT_SOURCE_DIR})
# AddFormatTarget(model_test ${CMAKE_CURRENT_SOURCE_DIR})
# AddFormatTarget(observation_test ${CMAKE_CURRENT_SOURCE_DIR})
# AddFormatTarget(obsoperator_test ${CMAKE_CURRENT_SOURCE_DIR})
#AddFormatTarget(geometry_test ${CMAKE_CURRENT_SOURCE_DIR})