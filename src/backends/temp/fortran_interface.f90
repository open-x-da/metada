module fortran_interface
  use iso_c_binding
  use mod_misc
  use mod_csp
  use mod_csp_init
  USE mod_mpi_interfaces
  USE mod_mpi_variables
  USE mod_mpi_csp_io
  USE mod_mpi_test
  !YY
  USE mitice
  USE mitice_parameters
  USE mitice_utility
  USE mitice_init
  USE mitice_ave
  use mod_csp_basic
  implicit none
  
  ! Add module variables for logging configuration
  logical, parameter :: log_enabled = .true.
  logical, parameter :: log_mpi_rank_only_0 = .true.
  character(len=24), parameter :: log_prefix = "[fortran_interface] "
    
  contains
  
  ! =================================================================
  ! Utility Functions
  ! =================================================================
  
  ! Central logging function for consistent output formatting
  subroutine log_message(message, force_all_ranks)
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: force_all_ranks
    logical :: print_message
    character(len=512) :: full_message
    
    ! Only proceed if logging is enabled
    if (.not. log_enabled) return
    
    ! Determine whether to print based on rank
    if (present(force_all_ranks)) then
      print_message = force_all_ranks .or. .not. log_mpi_rank_only_0 .or. mpi_rank == 0
    else
      print_message = .not. log_mpi_rank_only_0 .or. mpi_rank == 0
    end if
    
    if (print_message) then
      ! Format the message with prefix and rank info
      if (log_mpi_rank_only_0 .and. mpi_rank == 0) then
        write(full_message, '(A,A)') log_prefix, trim(message)
      else
        write(full_message, '(A,A,I0,A,A)') log_prefix, "Rank ", mpi_rank, ": ", trim(message)
      end if
      
      ! Print the message
      write(*, *) trim(full_message)
    end if
  end subroutine log_message
  
  ! Function to log debug messages (can be disabled separately if needed)
  subroutine log_debug(message, force_all_ranks)
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: force_all_ranks
    character(len=256) :: debug_message
    
    call log_message(message, force_all_ranks)
  end subroutine log_debug
  
  ! Function to log error messages (always displayed)
  subroutine log_error(message)
    character(len=*), intent(in) :: message
    character(len=256) :: error_message
    
    write(error_message, '(A,A)') "ERROR: ", trim(message)
    call log_message(error_message, .true.)
  end subroutine log_error
  
  ! =================================================================
  ! Interface Functions
  ! =================================================================
  
  subroutine c_mpi_process_init(comm) bind(C, name="c_mpi_process_init")
    ! Initialize MPI process, set communicator and call original initialization function
    integer(c_int), value :: comm
    character(len=256) :: msg
    
    ! Set MPI communicator
    mpi_comp_comm = comm
    
    ! Call original mpi_process_init function
    call mpi_process_init
    
    ! Add some debug output
    write(msg, '(A,I0)') "After mpi_process_init, mpi_io_procs=", mpi_io_procs
    call log_message(msg)

  end subroutine c_mpi_process_init
  
  subroutine c_read_namelist() bind(C, name="c_read_namelist")
    ! Read model configuration file
    call misc_namelist_read()
  end subroutine

  ! --- Sea Ice Module Related ---
  subroutine c_init_mitice_run_info() bind(C, name="c_init_mitice_run_info")
    ! Initialize sea ice module run information
    call mitice_run_info_open()
  end subroutine
  
  subroutine c_read_mitice_params() bind(C, name="c_read_mitice_params")
    ! Read sea ice module parameters
    call mitice_read_params()
  end subroutine

  subroutine c_allocate_mitice() bind(C, name="c_allocate_mitice")
    ! Allocate memory for sea ice module
    call mitice_init_allocate()
  end subroutine
  
  subroutine c_init_mitice_fixed() bind(C, name="c_init_mitice_fixed")
    ! Initialize sea ice module fixed parameters
    call mitice_init_fixed()
  end subroutine
  
  subroutine c_update_gpu_seaice_parameters() bind(C, name="c_update_gpu_seaice_parameters")
    ! Update GPU sea ice parameters
    call gpu_seaice_parameters_update()
  end subroutine
  
  subroutine c_init_mitice_vars() bind(C, name="c_init_mitice_vars")
    ! Initialize sea ice module variables
    call mitice_init_vars()
  end subroutine
  
  ! --- CSP Computation Related ---
  subroutine c_init_csp() bind(C, name="c_init_csp")
    ! Initialize CSP model
    call csp_init()
  end subroutine
  
  subroutine c_init_csp_asm() bind(C, name="c_init_csp_asm")
    ! Initialize CSP assimilation module
    call csp_init_asm()
  end subroutine
  
  subroutine c_init_gpu() bind(C, name="c_init_gpu")
    ! Initialize GPU variables
    call GPU_VARS_INIT()
  end subroutine
  
  subroutine c_call_csp(c_iter) bind(C, name="c_call_csp")
    ! Execute single step CSP computation
    integer(c_int), value, intent(in) :: c_iter
    character(len=256) :: msg

    myIter = c_iter
    
    ! Execute CSP computation
    write(msg, '(A,I0)') "myIter:", myIter
    call log_message(msg)

    call csp()
  end subroutine
  
  subroutine c_call_csp2() bind(C, name="c_call_csp2")
    ! Execute CSP main loop computation
    ! Execute loop from nIter0 to nIterMax
    do myIter = nIter0, nIterMax
        IF (mpi_rank == 0) write(*,*) "myIter:", myIter
        call csp()
    end do
  end subroutine


  subroutine c_finalize_mitice() bind(C, name="c_finalize_mitice")
    ! Clean up sea ice module resources
    call mitice_run_info_close()
    if (SEAICE_taveFreq > 0.0_wp) then
      call mitice_ave_release()
    endif
  end subroutine
  
  subroutine c_finalize_mpi() bind(C, name="c_finalize_mpi")
    ! Clean up MPI resources
    call mpi_final_operations()
  end subroutine

  ! ==================================================================
  ! II. Get or Send MPI Rank and Configuration Flags  
  ! ==================================================================
  
  subroutine c_get_mpi_rank(rank) bind(C, name="c_get_mpi_rank")
    ! Get current process MPI rank
    integer(c_int), intent(out) :: rank
    
    ! Return current process rank
    rank = mpi_rank
  end subroutine c_get_mpi_rank

  subroutine c_get_config_flags(mitice, restart, assim, init_iter, max_iter) bind(C, name="c_get_config_flags")
    ! Get model configuration flags and parameters
    logical(c_bool), intent(out) :: mitice, restart, assim
    integer(c_int), intent(out) :: init_iter, max_iter
    
    ! Get all configuration flags and parameters at once
    mitice = mitice_on
    restart = restart_in
    assim = assim_in
    init_iter = nIter0
    max_iter = nIterMax
  end subroutine
  
  function c_is_mitice_enabled() bind(C, name="c_is_mitice_enabled") result(enabled)
    ! Check if sea ice module is enabled
    logical(c_bool) :: enabled
    enabled = mitice_on
  end function
  
  function c_is_restart_enabled() bind(C, name="c_is_restart_enabled") result(enabled)
    ! Check if restart feature is enabled
    logical(c_bool) :: enabled
    enabled = restart_in
  end function
  
  function c_is_assim_enabled() bind(C, name="c_is_assim_enabled") result(enabled)
    ! Check if assimilation feature is enabled
    logical(c_bool) :: enabled
    enabled = assim_in
  end function
  
  function c_is_compute_process() bind(C, name="c_is_compute_process") result(is_comp)
    ! Check if current process is a compute process
    logical(c_bool) :: is_comp
    is_comp = mpi_rank < mpi_comp_procs
  end function  
  
  ! ==================================================================
  ! III. Macom Model Data and IO Management Related Functions
  ! ==================================================================
  
  ! --- Information Exchange Related ---
  subroutine c_send_info_to_io() bind(C, name="c_send_info_to_io")
    ! Send initial information to IO processes
    if (mpi_rank == 0) then
      call mpi_send_info_comp_to_io()
    end if
  end subroutine
  
  subroutine c_send_second_info_to_io() bind(C, name="c_send_second_info_to_io")
    ! Send second batch of information to IO processes
    if (mpi_rank == 0) then
      call mpi_send_2nd_info_comp_to_io()
    end if
  end subroutine
  
  ! --- IO Process Management ---
  subroutine c_init_io() bind(C, name="c_init_io")
    ! Initialize IO processes
    call mpi_csp_io_init()
  end subroutine
  
  subroutine c_csp_io_main() bind(c, name="c_csp_io_main")
    ! Run IO process main loop
    call mpi_csp_io_main()
  end subroutine
  
  subroutine c_send_output() bind(C, name="c_send_output")
    ! Send output data
    call mpi_csp_io_send_output()
  end subroutine
  
  subroutine c_read_restart() bind(C, name="c_read_restart")
    ! Read restart files
    call mpi_csp_io_restart_read()
  end subroutine
  
  ! --- Run Information Management ---
  subroutine c_open_misc_run_info() bind(C, name="c_open_misc_run_info")
    ! Open run information file
    call misc_run_info_open()
  end subroutine
  
  subroutine c_close_misc_run_info() bind(C, name="c_close_misc_run_info")
    ! Close run information file
    call misc_run_info_close()
  end subroutine
  
  ! ==================================================================
  ! IV. Fortran and C++ Data exchange related functions in Grid Interface
  ! ==================================================================

  ! --- Grid Data Access ---
  subroutine c_get_grid_dimensions(verticalLevels, verticalLevelsP1, dimensions, &
    totalPoints, computePoints, totalVortPoints, &
    computeVortPoints, boundaryPoints) &
    bind(C, name="c_get_grid_dimensions")
    ! Get grid dimension information
    integer(c_int), intent(out) :: verticalLevels, verticalLevelsP1, dimensions
    integer(c_int), intent(out) :: totalPoints, computePoints, totalVortPoints
    integer(c_int), intent(out) :: computeVortPoints, boundaryPoints

    ! Assign from module variables to interface variables
    verticalLevels = nk
    verticalLevelsP1 = nkp1
    dimensions = ni
    totalPoints = nl
    computePoints = nlpb
    totalVortPoints = nlz
    computeVortPoints = nlpbz
    boundaryPoints = nlbdy
  end subroutine

  subroutine c_get_vertical_grid_vars(centerDepthP_ptr, surfaceDepthP_ptr, &
                                   centerThicknessP_ptr, surfaceThicknessP_ptr, &
                                   centerDepthM_ptr, surfaceDepthM_ptr, &
                                   levelDensity_ptr) bind(C, name="c_get_vertical_grid_vars")
    ! Get vertical grid variables
    type(c_ptr), intent(out) :: centerDepthP_ptr, surfaceDepthP_ptr
    type(c_ptr), intent(out) :: centerThicknessP_ptr, surfaceThicknessP_ptr
    type(c_ptr), intent(out) :: centerDepthM_ptr, surfaceDepthM_ptr
    type(c_ptr), intent(out) :: levelDensity_ptr
    
    ! Directly pass array pointers
    centerDepthP_ptr = c_loc(rC)
    surfaceDepthP_ptr = c_loc(rF)
    centerThicknessP_ptr = c_loc(drC)
    surfaceThicknessP_ptr = c_loc(drF)
    centerDepthM_ptr = c_loc(rC_z)
    surfaceDepthM_ptr = c_loc(rF_z)
    levelDensity_ptr = c_loc(rhoLev)
  end subroutine
  
  subroutine c_get_horizontal_grid_vars(tracerDx_ptr, tracerDy_ptr, &
                                     uPointDx_ptr, uPointDy_ptr, &
                                     vPointDx_ptr, vPointDy_ptr, &
                                     vortexDx_ptr, vortexDy_ptr) bind(c, name="c_get_horizontal_grid_vars")
    ! Get horizontal grid variables
    type(c_ptr), intent(out) :: tracerDx_ptr, tracerDy_ptr
    type(c_ptr), intent(out) :: uPointDx_ptr, uPointDy_ptr
    type(c_ptr), intent(out) :: vPointDx_ptr, vPointDy_ptr
    type(c_ptr), intent(out) :: vortexDx_ptr, vortexDy_ptr
    
    ! Directly pass array pointers
    tracerDx_ptr = c_loc(dxC)
    tracerDy_ptr = c_loc(dyC)
    uPointDx_ptr = c_loc(dxW)
    uPointDy_ptr = c_loc(dyW)
    vPointDx_ptr = c_loc(dxS)
    vPointDy_ptr = c_loc(dyS)
    vortexDx_ptr = c_loc(dxZ)
    vortexDy_ptr = c_loc(dyZ)
  end subroutine

  subroutine c_get_geographic_vars(tracerLat_ptr, tracerLon_ptr, &
                                  uPointLat_ptr, uPointLon_ptr, &
                                  vPointLat_ptr, vPointLon_ptr, &
                                  coriolisParam_ptr) bind(c, name="c_get_geographic_vars")
    ! Get geographic variable information
    type(c_ptr), intent(out) :: tracerLat_ptr, tracerLon_ptr
    type(c_ptr), intent(out), optional :: uPointLat_ptr, uPointLon_ptr
    type(c_ptr), intent(out), optional :: vPointLat_ptr, vPointLon_ptr
    type(c_ptr), intent(out), optional :: coriolisParam_ptr
    
    ! Directly pass array pointers
    tracerLat_ptr = c_loc(latC)
    tracerLon_ptr = c_loc(lonC)
  end subroutine
  
  subroutine c_get_grid_mask_vars(tracerMask_ptr, uPointMask_ptr, vPointMask_ptr, vortexMask_ptr, &
                               tracerFactor_ptr, uPointFactor_ptr, vPointFactor_ptr) bind(c, name="c_get_grid_mask_vars")
    ! Get grid mask variables
    type(c_ptr), intent(out) :: tracerMask_ptr, uPointMask_ptr, vPointMask_ptr, vortexMask_ptr
    type(c_ptr), intent(out) :: tracerFactor_ptr, uPointFactor_ptr, vPointFactor_ptr
    
    ! Pass mask and factor pointers
    tracerMask_ptr = c_loc(maskC(1,1))
    uPointMask_ptr = c_loc(maskW(1,1))
    vPointMask_ptr = c_loc(maskS(1,1))
    vortexMask_ptr = c_loc(maskZ(1,1))
    tracerFactor_ptr = c_loc(h0FacC(1,1))
    uPointFactor_ptr = c_loc(h0FacW(1,1))
    vPointFactor_ptr = c_loc(h0FacS(1,1))
  end subroutine

  subroutine c_get_grid_area_vars(tracerArea_ptr, uPointArea_ptr, vPointArea_ptr, vortexArea_ptr, &
                                 metricQ11_ptr, metricQ12_ptr, metricQ22_ptr) bind(c, name="c_get_grid_area_vars")
    ! Get grid area variables
    type(c_ptr), intent(out) :: tracerArea_ptr, uPointArea_ptr, vPointArea_ptr, vortexArea_ptr
    type(c_ptr), intent(out), optional :: metricQ11_ptr, metricQ12_ptr, metricQ22_ptr
    
    ! Directly pass array pointers
    tracerArea_ptr = c_loc(rAc)
    uPointArea_ptr = c_loc(rAw)
    vPointArea_ptr = c_loc(rAs)
    vortexArea_ptr = c_loc(rAz)
  end subroutine
  
  subroutine c_get_velocity_fields(uField_ptr, vField_ptr, wField_ptr) bind(c, name="c_get_velocity_fields")
    ! Get velocity field data
    type(c_ptr), intent(out) :: uField_ptr, vField_ptr, wField_ptr
    character(len=256) :: msg
    integer :: dim1, dim2

    ! Check if arrays are allocated
    if (allocated(uFld) .and. allocated(vFld) .and. allocated(wFld)) then
        ! Get dimensions properly
        
        ! ! Store dimensions in local variables first
        ! dim1 = size(uFld, 1)
        ! dim2 = size(uFld, 2)
        
        ! ! Then use the variables in write statement
        ! write(msg, '(A,I0,A,I0,A)') "Velocity fields dimensions: ", &
        !       dim1, " x ", dim2, " (nlpb x nk)"
        ! call log_message(msg)
        
        ! Use address of first element of arrays
        uField_ptr = c_loc(uFld(1,1))
        vField_ptr = c_loc(vFld(1,1))
        wField_ptr = c_loc(wFld(1,1))
        
        ! ! Debug: Print first value from each velocity field
        ! if (size(uFld) > 0) then
        !     write(msg, '(A,E14.6,A,E14.6,A,E14.6)') "First velocity values: u=", &
        !           uFld(1,1), ", v=", vFld(1,1), ", w=", wFld(1,1)
        !     call log_message(msg)
        !     write(msg, '(A,E14.6,A,E14.6,A,E14.6)') "First velocity values: u=", &
        !           uFld(1,2), ", v=", vFld(1,2), ", w=", wFld(1,2)
        !     call log_message(msg)
        ! end if
        
        ! call log_message("Velocity field pointers passed to C++")
    else
        ! ! Return null pointers if arrays not allocated
        ! call log_message("ERROR: Velocity field arrays not allocated!")
        uField_ptr = c_null_ptr
        vField_ptr = c_null_ptr
        wField_ptr = c_null_ptr
    end if
  end subroutine
  
  subroutine c_get_temperature_salinity(temperatureField_ptr, salinityField_ptr) bind(c, name="c_get_temperature_salinity")
    ! Get temperature and salinity field data
    type(c_ptr), intent(out) :: temperatureField_ptr, salinityField_ptr
    character(len=256) :: msg
    
    ! Check if arrays are allocated
    if (allocated(tFld)) then
      ! Get pointers to temperature and salinity fields
      ! tFld(1,1,1) is temperature, tFld(1,1,2) is salinity
      temperatureField_ptr = c_loc(tFld(1,1,1))
      salinityField_ptr = c_loc(tFld(1,1,2))
      
      ! ! Debug: Print first value from temperature and salinity fields
      ! if (size(tFld) > 0) then
      !   write(msg, '(A,E14.6,A,E14.6)') "First temperature and salinity values: T=", &
      !         tFld(1,1,1), ", S=", tFld(1,1,2)
      !   call log_message(msg)
        
      !   if (size(tFld, 2) > 1) then
      !     write(msg, '(A,E14.6,A,E14.6)') "Second layer values: T=", &
      !           tFld(1,2,1), ", S=", tFld(1,2,2)
      !     call log_message(msg)
      !   end if
      ! end if
      
      ! call log_message("Temperature and salinity field pointers passed to C++")
    else
      ! ! Return null pointers if arrays not allocated
      ! call log_message("ERROR: Temperature and salinity field arrays not allocated!")
      temperatureField_ptr = c_null_ptr
      salinityField_ptr = c_null_ptr
    end if
  end subroutine


  subroutine c_set_velocity_fields(uField, vField, wField, size, i_start, j_start, ni, nj) bind(c, name="c_set_velocity_fields")
    ! Set velocity field data with options for global or local modification
    ! Parameters:
    ! - uField, vField, wField: Velocity components (required)
    ! - size: Total size of arrays (required for validation)
    ! - i_start, j_start, ni, nj: Optional region parameters
    !   If not provided or <= 0, defaults to global modification
    real(c_double), intent(in) :: uField(*), vField(*), wField(*)
    integer(c_int), value, intent(in) :: size
    integer(c_int), intent(in), optional :: i_start, j_start, ni, nj
    character(len=256) :: msg
    integer :: i, j, idx, nlp, nk_local
    integer :: i_begin, j_begin, i_end, j_end, count
    logical :: global_mode
    
    ! Check if arrays are allocated
    if (allocated(uFld) .and. allocated(vFld) .and. allocated(wFld)) then
      ! Calculate dimensions
      nlp = nlpb  ! Number of horizontal points
      nk_local = nk  ! Number of vertical levels
      
      ! Check if size matches our expectations
      if (size /= nlp * nk_local) then
        write(msg, '(A,I0,A,I0,A,I0,A,I0)') "ERROR: Size mismatch! Got ", size, &
              " elements, expected ", nlp * nk_local, " (", nlp, " × ", nk_local, ")"
        call log_message(msg, .true.)
        return
      end if
      
      ! Set default region (global)
      i_begin = 1
      j_begin = 1
      i_end = nlp
      j_end = nk_local
      
      ! Determine if we're in global or local modification mode based on parameters
      global_mode = .true.
      
      ! Override with provided values if present and valid
      if (present(i_start) .and. i_start > 0) then
        i_begin = i_start
        global_mode = .false.
      end if
      
      if (present(j_start) .and. j_start > 0) then
        j_begin = j_start
        global_mode = .false.
      end if
      
      ! Handle region extents
      if (present(ni) .and. ni > 0) then
        i_end = min(i_begin + ni - 1, nlp)
      end if
      
      if (present(nj) .and. nj > 0) then
        j_end = min(j_begin + nj - 1, nk_local)
      end if
      
      ! Log what we're doing
      if (global_mode) then
        call log_message("Copying ALL velocity data from C++ to Fortran arrays")
      else
        write(msg, '(A,I0,A,I0,A,I0,A,I0)') "Modifying LOCAL region: i=", i_begin, ":", i_end, &
                                           ", j=", j_begin, ":", j_end
        call log_message(msg)
      end if
      
      ! Update the data - same loop for both global and local cases,
      ! just with different bounds
      count = 0
      do j = j_begin, j_end
        do i = i_begin, i_end
          ! Calculate 0-based C++ index from our region
          idx = (j-1) * nlp + (i-1) + 1
          
          ! Safety check for index
          if (idx > 0 .and. idx <= size) then
            uFld(i,j) = uField(idx)
            vFld(i,j) = vField(idx)
            wFld(i,j) = wField(idx)
            count = count + 1
          end if
        end do
      end do
      
      ! Log first few values
      write(msg, '(A,E14.6,A,E14.6,A,E14.6)') "First velocity values after update: u=", &
            uFld(1,1), ", v=", vFld(1,1), ", w=", wFld(1,1)
      call log_message(msg)
      
      write(msg, '(A,I0,A)') "Velocity field data updated successfully (", count, " points modified)"
      call log_message(msg)
    else
      call log_message("ERROR: Velocity field arrays not allocated, cannot update!", .true.)
    end if
  end subroutine

end module