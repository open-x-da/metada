MODULE mod_mpi_test_variables 
  use mod_misc_basic
  IMPLICIT NONE
   ! for mpi test function
   INTEGER, PUBLIC, PARAMETER :: MPI_TEST_UNIT = 93   ! unit for outputing performance test results for computing processes
   INTEGER, PUBLIC, PARAMETER :: MPI_TEST_IO_UNIT = 94   ! unit for outputing performance test results for io processes
   CHARACTER(LEN = 10), PUBLIC :: mpi_time ! use time to name test file
   CHARACTER(LEN = lc), PUBLIC :: mpi_test_filename ! file name of test file
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_start  ! used for counting program runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_end  ! used for counting program runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_rate  ! used for counting program runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_count_max  ! used for counting program runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_compute_start  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_compute_end  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_init_start  ! used for counting initializing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_init_end  ! used for counting initializing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_communicate_start  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_communicate_end  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_start  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_end  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_read_start  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_read_end  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_write_start  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_write_end  ! used for counting computing runtime with SYSTEM_CLOCK function
   
   ! for mpi debug function
   INTEGER, PUBLIC :: mpi_debug_size  ! test for size of integer and real in mpi
   INTEGER, PUBLIC :: mpi_int_type ! test for size of integer in fortran
   REAL(KIND = dp) :: mpi_real_type ! test for size of real in fortran
   INTEGER, PUBLIC, PARAMETER :: MPI_DEBUG_UNIT = 92   ! unit for outputing to debug file
   CHARACTER(LEN = lc), PUBLIC :: mpi_debug_filename ! file name of debug file
   CHARACTER(LEN = lc), PUBLIC :: mpi_var_filename ! file name of debug file

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :, :) :: mpi_variable_3d_nk_2 ! test for size of real in fortran
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_variable_2d_nk ! test for size of real in fortran

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_variable_1d ! test for size of real in fortran
   REAL(dp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_dp_1d ! test for size of real in fortran
   REAL(dp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_dp_2d_nk ! test for size of real in fortran

END MODULE

