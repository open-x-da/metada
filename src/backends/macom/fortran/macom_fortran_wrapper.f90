module macom_fortran_wrapper
  use iso_c_binding
  use macom_logger
  ! Assuming your actual model logic is in other modules, e.g., mod_csp, mod_misc etc.
  ! You will need to add 'use' statements for any modules that these wrapper subroutines call.
  ! For example:
  use mod_misc
  use mod_csp
  use mod_csp_init
  use mod_mpi_interfaces
  use mod_mpi_variables
  use mod_mpi_csp_io
  use mod_mpi_test
  ! Conditionally include mitice modules based on compilation flags

#ifdef MACOM_REGION_MODE
  use mod_macom_final
  use mod_csp_basic
#endif

#ifdef SEAICE_ITD
  use mitice
  use mitice_parameters
  use mitice_utility
  use mitice_init
  use mitice_ave
#endif
  implicit none

  private ! Default to private, only expose BIND(C) interfaces

  ! ! Publicly expose only the C-bindable procedures
  ! public :: c_macom_initialize_mpi
  ! public :: c_macom_get_mpi_rank
  ! public :: c_macom_get_mpi_size
  ! public :: c_macom_read_namelist
  ! public :: c_macom_finalize_mpi

  ! Example: Variables from a shared module that might be needed by the wrappers.
  ! Ensure these are properly 'use'd from their respective modules.
  ! integer :: mpi_comp_comm
  ! integer :: mpi_rank
  ! logical :: mitice_on
  ! logical :: restart_in
  ! logical :: assim_in
  ! integer :: nIter0, nIterMax
  ! integer :: myIter ! Current iteration step

contains

  !-----------------------------------------------------------------------------
  ! I. MPI Management
  !-----------------------------------------------------------------------------
  subroutine c_macom_initialize_mpi(comm_cpp) &
    bind(C, name="c_macom_initialize_mpi")
    integer(C_INT), value, intent(in) :: comm_cpp
    integer :: ierr

    ! Convert Cpp communicator to Fortran communicator
    ! wrapper_mpi_comm = comm_cpp

    ! Set MPI communicator
    ! mpi_comp_comm = comm_cpp

#ifdef MACOM_REGION_MODE
    ! Initialize loop variables for region mode
    macom_loop_i = 1
    macom_loop_end = 1
#endif

    ! Call MACOM's MPI initialization
    ! call init_mpi()
    call mpi_process_init

  end subroutine c_macom_initialize_mpi

  subroutine c_macom_get_mpi_rank(rank) &
    bind(C, name="c_macom_get_mpi_rank")
    integer(C_INT), intent(out) :: rank
    ! Return current process rank
    rank = mpi_rank
    ! rank = wrapper_mpi_rank
  end subroutine c_macom_get_mpi_rank

  subroutine c_macom_get_mpi_size(size) &
    bind(C, name="c_macom_get_mpi_size")
    integer(C_INT), intent(out) :: size
    ! Return current process size
    size = mpi_procs
    ! size = wrapper_mpi_size
  end subroutine c_macom_get_mpi_size

  subroutine c_macom_finalize_mpi() &
    bind(C, name="c_macom_finalize_mpi")
    ! Use conditional compilation to handle different finalization methods
#ifdef MACOM_GLOBAL_MODE
    ! Global mode: use mpi_final_operations
    call mpi_final_operations()
    if (mpi_rank .eq. 0) then
      call macom_log_info("FortranWrapper", &
                "c_macom_finalize_mpi: used mpi_final_operations (global mode)")
    end if
#endif

#ifdef MACOM_REGION_MODE
    ! Region mode: use region-specific finalization
    ! call macom_final_mpi
    call mpi_netcdf_read_finalize
    call MPI_FINALIZE(mpi_err)
    if (mpi_rank .eq. 0) then
      call macom_log_info("FortranWrapper", &
                          "c_macom_finalize_mpi: used region mode finalization")
    end if
#endif

  end subroutine c_macom_finalize_mpi

  !-----------------------------------------------------------------------------
  ! II. Model Configuration and Initialization
  !-----------------------------------------------------------------------------
  subroutine c_macom_read_namelist() &
    bind(C, name="c_macom_read_namelist")
    ! Reads the model configuration (namelist).
    call misc_namelist_read()
    if (mpi_rank .eq. 0) then
      call macom_log_info("FortranWrapper", "c_macom_read_namelist called")
    end if
  end subroutine c_macom_read_namelist

  subroutine c_get_macom_config_flags(mitice, restart, assim, &
                   init_iter, max_iter) bind(c, name='c_get_macom_config_flags')
    ! Get all configuration flags and parameters at once
    logical(C_BOOL), intent(out) :: mitice, restart, assim
    integer(C_INT), intent(out) :: init_iter, max_iter

    ! Use actual mitice_on value from namelist
#ifdef SEAICE_ITD
    mitice = mitice_on
#else
    mitice = .false.
#endif

    restart = restart_in
    assim = assim_in
    init_iter = nIter0
    max_iter = nIterMax

    if (mpi_rank .eq. 0) then
      call macom_log_info("FortranWrapper", "c_get_macom_config_flags called")
    end if
  end subroutine c_get_macom_config_flags

  subroutine c_macom_initialize_mitice() &
    bind(C, name="c_macom_initialize_mitice")
#ifdef SEAICE_ITD
    if (mitice_on) then
      ! Open runtime file for screen output and model parameter summary
      call mitice_run_info_open()
      ! Read seaice parameter namelist /namelist.mitice/
      call mitice_read_params()
      if (mpi_rank .eq. 0) then
        call macom_log_info("FortranWrapper", &
                            "c_macom_initialize_mitice called")
      end if
    end if
#else
    ! For builds without sea ice support, this is a no-op
    if (mpi_rank .eq. 0) then
      call macom_log_info("FortranWrapper", &
              "c_macom_initialize_mitice called (no-op - sea ice not compiled)")
    end if
#endif
  end subroutine c_macom_initialize_mitice

  !-----------------------------------------------------------------------------
  ! Custom: Expose mpi_send_info_comp_to_io and misc_run_info_open to C++
  !-----------------------------------------------------------------------------
  subroutine c_macom_mpi_send_info_comp_to_io() &
    bind(C, name="c_macom_mpi_send_info_comp_to_io")
    if (mpi_rank .eq. 0) then
      call mpi_send_info_comp_to_io()
    end if
  end subroutine c_macom_mpi_send_info_comp_to_io

  subroutine c_macom_misc_run_info_open() &
    bind(C, name="c_macom_misc_run_info_open")
    call misc_run_info_open()
  end subroutine c_macom_misc_run_info_open

  !-----------------------------------------------------------------------------
  ! CSP Initialization (C++ interface)
  !-----------------------------------------------------------------------------
  subroutine c_macom_init_csp() &
    bind(C, name="c_macom_init_csp")
    ! send information to io processes from computational processes
    if (mpi_rank .eq. 0) then
      call macom_log_info("FortranWrapper", "c_macom_init_csp")
    end if
    call csp_init()
  end subroutine c_macom_init_csp

  !-----------------------------------------------------------------------------
  ! Mitice full initialization (C++ interface)
  !-----------------------------------------------------------------------------
  subroutine c_macom_mitice_init_all() &
    bind(C, name="c_macom_mitice_init_all")
#ifdef SEAICE_ITD
    if (mitice_on) then
      call mitice_init_allocate !allocate host and device variables
      call mitice_init_fixed !prepare a few parameters
      call gpu_seaice_parameters_update !upload parameters to GPU device
      call mitice_init_vars !initialized model variables and read initial fields
    end if
#else
    ! For builds without sea ice support, this is a no-op
    if (mpi_rank .eq. 0) then
      call macom_log_info("FortranWrapper", &
                "c_macom_mitice_init_all called (no-op - sea ice not compiled)")
    end if
#endif
  end subroutine c_macom_mitice_init_all

  !-----------------------------------------------------------------------------
  ! Restart and assimilation (C++ interface)
  !-----------------------------------------------------------------------------
  subroutine c_macom_restart_and_assim() &
    bind(C, name="c_macom_restart_and_assim")
    if (restart_in) then
      call mpi_csp_io_restart_read
      if (assim_in) then
        call csp_init_asm
      end if
    end if

    if (mpi_rank .eq. 0) then
      ! send additional information from comp communicator to io communicator
      call mpi_send_2nd_info_comp_to_io
    end if

    ! ---- OPENACC: UPDATE VARS TO GPU FOR TIME INTEGRATION
    call GPU_VARS_INIT
    ! ----

    if (.not. restart_in) then
      ! send output to IO process
      call mpi_csp_io_send_output
    end if

  end subroutine c_macom_restart_and_assim

  !-----------------------------------------------------------------------------
  ! CSP single step (C++ interface)
  !-----------------------------------------------------------------------------
  subroutine c_macom_run_csp_step() &
    bind(C, name="c_macom_run_csp_step")
    do myIter = nIter0, nIterMax

      if (mpi_rank .eq. 0) write (*, *) "myIter:", myIter

      call csp

    end do

    call misc_run_info_close
  end subroutine c_macom_run_csp_step

  !-----------------------------------------------------------------------------
  ! IO process main (C++ interface)
  !-----------------------------------------------------------------------------
  subroutine c_macom_csp_io_main() &
    bind(C, name="c_macom_csp_io_main")
    call mpi_csp_io_init
    call mpi_csp_io_main
  end subroutine c_macom_csp_io_main

  !-----------------------------------------------------------------------------
  ! Mitice finalize (C++ interface) - Original function for global mode
  !-----------------------------------------------------------------------------
  subroutine c_macom_finalize_mitice() &
    bind(C, name="c_macom_finalize_mitice")
#ifdef SEAICE_ITD
    if (mitice_on) then
      call mitice_run_info_close !close seaice running message file
      if (SEAICE_taveFreq .gt. 0.0) then
        call mitice_ave_release
      end if
    end if
#else
    ! For builds without sea ice support, this is a no-op
    if (mpi_rank .eq. 0) then
      call macom_log_info("FortranWrapper", &
                "c_macom_finalize_mitice called (no-op - sea ice not compiled)")
    end if
#endif
  end subroutine c_macom_finalize_mitice

  !-----------------------------------------------------------------------------
  ! Utility functions
  !-----------------------------------------------------------------------------
  function c_int_to_string(i) result(str)
    integer(C_INT), intent(in) :: i
    character(len=20) :: str

    write (str, '(I0)') i
  end function c_int_to_string

end module macom_fortran_wrapper
