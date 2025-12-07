setup_test_executable(datetime_test DateTimeTest.cpp)
add_test_format_target(datetime_test ${CMAKE_CURRENT_SOURCE_DIR})

setup_test_executable(duration_test DurationTest.cpp)
add_test_format_target(duration_test ${CMAKE_CURRENT_SOURCE_DIR})

# Setup MPI context test with MPI support
# This will automatically wrap tests with mpirun when ctest is run
if(METADA_USE_MPI AND MPI_FOUND)
  # Create the executable first
  add_executable(mpicontext_test MPIContextTest.cpp)
  
  target_link_libraries(mpicontext_test
    PRIVATE
      metada::traits
      metada::framework::adapters
      metada::framework::runs
      metada::base
      metada::framework::algorithms
      metada::backends::gmock
  )
  
  # Use MPI-aware test setup (this adds coverage and discovers tests)
  # CROSSCOMPILING_EMULATOR will automatically wrap all test commands with mpirun
  metada_add_test_with_coverage_mpi(mpicontext_test 4)
  add_test_format_target(mpicontext_test ${CMAKE_CURRENT_SOURCE_DIR})
else()
  # Fall back to regular test setup if MPI is not available
  setup_test_executable(mpicontext_test MPIContextTest.cpp)
  add_test_format_target(mpicontext_test ${CMAKE_CURRENT_SOURCE_DIR})
endif()