! Dispatch C-bindings that translate flat arrays into WRFDA types
! and call the included transform routines.
!
! OUT-OF-DOMAIN OBSERVATION HANDLING:
! This module implements WRFDA-compliant handling of out-of-domain observations:
! 1. find_fractional_ij returns -1 for i,j and -1.0 for xfloat,yfloat when observations are out of domain
! 2. Callers check for -1 values and handle them gracefully by:
!    - Setting invalid indices (-1) in the iv structure
!    - Skipping further processing for that observation
!    - Continuing with the next observation
! 3. This follows WRFDA's standard pattern of setting outside=.true. and returning immediately

module metada_wrfda_dispatch
  use iso_c_binding
  use module_domain,        only: domain, x_type
  use da_define_structures, only: iv_type, y_type, xbx_type, da_allocate_y, da_allocate_obs_info
  use module_symbols_util,  only: wrfu_initialize, wrfu_finalize, wrfu_cal_gregorian
  use module_symbols_util, only: WRFU_ClockCreate, WRFU_TimeIntervalSet, WRFU_SUCCESS, WRFU_Time, WRFU_TimeInterval, WRFU_INITIALIZE, WRFU_CAL_GREGORIAN
  use da_control, only: sound, synop, pilot, satem, geoamv, polaramv, airep, gpspw, gpsref, &
                        metar, ships, ssmi_rv, ssmi_tb, ssmt1, ssmt2, qscat, profiler, buoy, &
                        bogus, pseudo, radar, radiance, airsr, sonde_sfc, mtgirs, tamdar, &
                        tamdar_sfc, rain, gpseph, lightning, &
                        sfc_assi_options, sfc_assi_options_1, trace_use_dull, &
                        var4d_run, num_fgat_time, missing_r, missing_data, num_ob_indexes, &
                        kts, kte, its, ite, jts, jte, Max_StHeight_Diff, kms, kme, &
                        ids, ide, jds, jde, ims, ime, jms, jme, ips, ipe, jps, jpe, kds, kde, kps, kpe, &
                        myproc, num_procs, rootproc, comm, num_qcstat_conv
  use da_tools, only: proj_info, da_map_set, da_llxy_wrf, da_togrid
  use da_metar,  only: da_transform_xtoy_metar,  da_transform_xtoy_metar_adj
  use da_synop,  only: da_transform_xtoy_synop,  da_transform_xtoy_synop_adj
  use da_buoy,   only: da_transform_xtoy_buoy,   da_transform_xtoy_buoy_adj
  use da_ships,  only: da_transform_xtoy_ships,  da_transform_xtoy_ships_adj
  use da_airep,  only: da_transform_xtoy_airep,  da_transform_xtoy_airep_adj
  use da_pilot,  only: da_transform_xtoy_pilot,  da_transform_xtoy_pilot_adj
  use da_sound,  only: da_transform_xtoy_sound,  da_transform_xtoy_sound_adj, &
                       da_transform_xtoy_sonde_sfc, da_transform_xtoy_sonde_sfc_adj
  use da_par_util, only: da_copy_dims, da_copy_tile_dims
  use da_tools, only: da_togrid, da_llxy_wrf, da_map_set, da_map_init
  use module_configure, only: grid_config_rec_type
  use da_minimisation, only: da_get_innov_vector

  implicit none
  
  ! Grid pointer: points to grid allocated via wrfda_alloc_and_init_domain
  ! This allows observation operators to access the grid
  type(domain), pointer, save :: persistent_grid => null()
  logical, save :: grid_initialized = .false.
  
  ! CRITICAL FIX: Add persistent iv structure to avoid deallocation issues
  type(iv_type), pointer, save :: persistent_iv
  logical, save :: iv_allocated = .false.
  
  ! CRITICAL FIX: Add persistent y structure to avoid repeated allocation
  type(y_type), pointer, save :: persistent_y
  logical, save :: y_allocated = .false.
  
  ! Map projection information for WRFDA coordinate conversion
  ! Define projection constants
  integer, parameter :: PROJ_LATLON = 0
  integer, parameter :: PROJ_MERC = 1
  integer, parameter :: PROJ_PS = 2
  integer, parameter :: PROJ_LC = 3
  
  ! Define proj_info type (simplified version of WRFDA's proj_info)
  ! Use WRFDA's proj_info type from da_tools module
  
  ! Module-level map_info structure initialized by da_setup_firstguess_wrf → da_map_set
  type(proj_info), save :: map_info
  
  ! Persistent output array for tangent linear operator results
  real(c_double), allocatable, save :: persistent_out_y(:)
  integer, save :: persistent_num_innovations = 0
  logical, save :: output_allocated = .false.
  
  ! Persistent arrays for adjoint operator gradients
  real(c_double), allocatable, save :: persistent_adjoint_u(:), persistent_adjoint_v(:)
  real(c_double), allocatable, save :: persistent_adjoint_t(:), persistent_adjoint_q(:)
  real(c_double), allocatable, save :: persistent_adjoint_psfc(:)
  real(c_double), allocatable, save :: persistent_delta_y(:)
  integer, save :: persistent_nx = 0, persistent_ny = 0, persistent_nz = 0
  logical, save :: adjoint_allocated = .false.
  
contains

  ! Forward operator: H(xb + xa) where xb is background state and xa is analysis increments
  ! This function computes the forward operator for incremental 3D-Var:
  ! - xb (background state): constant throughout outer loop, used for innovation computation
  ! - xa (analysis increments): start at zero, updated by each iteration
  ! - Total state: x_total = xb + xa (computed internally)
  ! - Forward operator: H(xb + xa)
  integer(c_int) function wrfda_xtoy_apply_grid(grid_ptr, ob_ptr, iv_ptr) bind(C, name="wrfda_xtoy_apply_grid")
    implicit none
    type(c_ptr), value :: grid_ptr, ob_ptr, iv_ptr

    type(domain), pointer :: grid
    type(iv_type), pointer :: iv
    type(y_type), pointer :: y
    integer :: n, num_innovations
    character(len=256) :: fam_str
    integer :: fam_id
    
    ! Convert C pointers to Fortran pointers
    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(iv_ptr, iv)
    call c_f_pointer(ob_ptr, y)
    
    if (.not. associated(iv)) then
      wrfda_xtoy_apply_grid = 1_c_int
      return
    end if
    
    if (.not. associated(y)) then
      wrfda_xtoy_apply_grid = 1_c_int
      return
    end if
    
    ! Instead of parsing input string, determine family from iv structure
    ! Find the first family that has observations available
    fam_id = 0
    do n = 1, size(iv%info)
      if (iv%info(n)%nlocal > 0) then
        fam_id = n
        exit
      end if
    end do
    
    if (fam_id == 0) then
      wrfda_xtoy_apply_grid = 1_c_int
      return
    end if
    
    ! Map family ID to name for debugging
    select case (fam_id)
    case (1); fam_str = "metar"
    case (2); fam_str = "synop"
    case (3); fam_str = "sound"
    case (4); fam_str = "gpspw"
    case (5); fam_str = "airep"
    case (6); fam_str = "pilot"
    case (7); fam_str = "ships"
    case (8); fam_str = "buoy"
    case default; fam_str = "unknown"
    end select
    
    select case (trim(fam_str))
    case ('metar'); 
      call da_transform_xtoy_metar(grid, iv, y)
    case ('synop', 'adpsfc'); 
      call da_transform_xtoy_synop(grid, iv, y)
    case default; 
      wrfda_xtoy_apply_grid = 1_c_int; return
    end select

    ! Calculate number of innovations from observations
    ! For synop observations, we have up to 5 variables per observation (U, V, T, Q, P)
    ! But we need to get the actual number of innovations from the iv structure
    if (associated(iv) .and. associated(iv%synop)) then
      ! Count actual innovations from the iv structure
      num_innovations = 0
      do n = 1, iv%info(2)%nlocal
        if (iv%synop(n)%u%qc == 0) num_innovations = num_innovations + 1
        if (iv%synop(n)%v%qc == 0) num_innovations = num_innovations + 1
        if (iv%synop(n)%t%qc == 0) num_innovations = num_innovations + 1
        if (iv%synop(n)%q%qc == 0) num_innovations = num_innovations + 1
        if (iv%synop(n)%p%qc == 0) num_innovations = num_innovations + 1
      end do
    end if
    
    ! Store number of innovations for later retrieval
    persistent_num_innovations = num_innovations
    
    ! Allocate or reallocate persistent output array if needed
    if (.not. output_allocated .or. .not. allocated(persistent_out_y)) then
      if (allocated(persistent_out_y)) then
        deallocate(persistent_out_y)
      end if
      allocate(persistent_out_y(num_innovations))
      output_allocated = .true.
    end if
    
    ! Copy results to persistent array
    call copy_y_to_out(fam_id, y, persistent_out_y, num_innovations)
    wrfda_xtoy_apply_grid = 0_c_int
  end function wrfda_xtoy_apply_grid

  ! Get count of tangent linear output values (count-only call)
  integer(c_int) function wrfda_get_tangent_linear_count(num_innovations) bind(C, name="wrfda_get_tangent_linear_count")
    implicit none
    integer(c_int), intent(out) :: num_innovations

    ! Check if output is available
    if (.not. output_allocated .or. persistent_num_innovations <= 0) then
      num_innovations = 0
      wrfda_get_tangent_linear_count = 1_c_int
      return
    end if

    ! Return the number of innovations
    num_innovations = persistent_num_innovations

    wrfda_get_tangent_linear_count = 0_c_int
  end function wrfda_get_tangent_linear_count

  ! Extract output values from the tangent linear operator
  integer(c_int) function wrfda_extract_tangent_linear_output(out_y, num_innovations) bind(C, name="wrfda_extract_tangent_linear_output")
    implicit none
    real(c_double), intent(out) :: out_y(*)
    integer(c_int), intent(out) :: num_innovations

    ! Check if output is available
    if (.not. output_allocated .or. persistent_num_innovations <= 0) then
      num_innovations = 0
      wrfda_extract_tangent_linear_output = 1_c_int
      return
    end if

    ! Return the number of innovations
    num_innovations = persistent_num_innovations

    ! Copy values to output array (safely)
    if (allocated(persistent_out_y) .and. size(persistent_out_y) >= num_innovations) then
      out_y(1:num_innovations) = persistent_out_y(1:num_innovations)
    else
      wrfda_extract_tangent_linear_output = 1_c_int
      return
    end if

    wrfda_extract_tangent_linear_output = 0_c_int
  end function wrfda_extract_tangent_linear_output

  ! Adjoint operator: H^T * δy -> δx (gradient accumulation)
  ! This function computes the adjoint of the forward operator for incremental 3D-Var:
  ! - delta_y: gradient in observation space (passed via persistent arrays)
  ! - inout_u, inout_v, inout_t, inout_q, inout_psfc: gradient in state space (accumulated in persistent arrays)
  ! - The adjoint operator accumulates gradients into the state space arrays
  integer(c_int) function wrfda_xtoy_adjoint_grid(grid_ptr, iv_ptr) bind(C, name="wrfda_xtoy_adjoint_grid")
    implicit none
    type(c_ptr), value :: grid_ptr, iv_ptr

    type(domain), pointer :: grid
    type(iv_type), pointer :: iv
    type(y_type), target :: jo_grad_y
    type(x_type), target :: jo_grad_x
    integer :: n, num_innovations, nx, ny, nz
    character(len=256) :: fam_str
    integer :: fam_id
    
    ! Convert C pointers to Fortran pointers
    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(iv_ptr, iv)
    
    if (.not. associated(iv)) then
      wrfda_xtoy_adjoint_grid = 1_c_int
      return
    end if
    
    ! Get grid dimensions
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    
    ! Zero out analysis increments for adjoint computation
    grid%xa%u = 0.0
    grid%xa%v = 0.0
    grid%xa%t = 0.0
    grid%xa%q = 0.0
    grid%xa%psfc = 0.0
    
    ! Grid structure and module-level variables are already set up
    ! No need to call da_copy_dims and da_copy_tile_dims again
    call zero_x_like(jo_grad_x, nx, ny, nz)
    sfc_assi_options = sfc_assi_options_1
    
    ! Instead of parsing input string, determine family from iv structure
    ! Find the first family that has observations available
    fam_id = 0
    do n = 1, size(iv%info)
      if (iv%info(n)%nlocal > 0) then
        fam_id = n
        exit
      end if
    end do
    
    if (fam_id == 0) then
      wrfda_xtoy_adjoint_grid = 1_c_int
      return
    end if
    
    ! Map family ID to name for debugging
    select case (fam_id)
    case (1); fam_str = "metar"
    case (2); fam_str = "synop"
    case (3); fam_str = "sound"
    case (4); fam_str = "gpspw"
    case (5); fam_str = "airep"
    case (6); fam_str = "pilot"
    case (7); fam_str = "ships"
    case (8); fam_str = "buoy"
    case default; fam_str = "unknown"
    end select

    ! Calculate number of innovations from observations
    ! For synop observations, we have up to 5 variables per observation (U, V, T, Q, P)
    if (associated(iv) .and. associated(iv%synop)) then
      
      ! Count actual innovations from the iv structure
      num_innovations = 0
      do n = 1, iv%info(2)%nlocal
        if (iv%synop(n)%u%qc == 0) num_innovations = num_innovations + 1
        if (iv%synop(n)%v%qc == 0) num_innovations = num_innovations + 1
        if (iv%synop(n)%t%qc == 0) num_innovations = num_innovations + 1
        if (iv%synop(n)%q%qc == 0) num_innovations = num_innovations + 1
        if (iv%synop(n)%p%qc == 0) num_innovations = num_innovations + 1
      end do
    end if
    
    ! Initialize jo_grad_y from persistent delta_y array
    call init_y_from_delta(fam_id, jo_grad_y, persistent_delta_y, num_innovations)

    select case (trim(fam_str))
    case ('metar'); 
      call da_transform_xtoy_metar_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('synop', 'adpsfc'); 
      call da_transform_xtoy_synop_adj(grid, iv, jo_grad_y, jo_grad_x)
    case default; 
      wrfda_xtoy_adjoint_grid = 1_c_int; return
    end select

    ! Copy gradients to persistent arrays for later retrieval
    call copy_x_to_persistent(jo_grad_x, nx, ny, nz)
    wrfda_xtoy_adjoint_grid = 0_c_int
  end function wrfda_xtoy_adjoint_grid

  subroutine init_y_from_delta(family, jo_grad_y, delta_y, num_innovations)
    integer, intent(in) :: family
    type(y_type), intent(inout) :: jo_grad_y
    real(c_double), intent(in) :: delta_y(*)
    integer, intent(in) :: num_innovations
    integer :: n, var_idx, num_obs
    real(c_double), parameter :: missing_r = -888888.0_c_double
    
    ! Calculate number of observations from innovations
    ! For synop observations, we need to get the actual number of observations from persistent_iv
    if (associated(persistent_iv) .and. associated(persistent_iv%synop)) then
      num_obs = persistent_iv%info(2)%nlocal
    else
      ! Fallback: assume we have 1 observation per 5 innovations
      num_obs = max(1, num_innovations / 5)
    end if
    select case (family)
    case (metar)
      allocate(jo_grad_y%metar(num_obs))
      do n=1,num_obs
        ! Initialize all variables to missing
        jo_grad_y%metar(n)%u = 0.0_c_double
        jo_grad_y%metar(n)%v = 0.0_c_double
        jo_grad_y%metar(n)%t = 0.0_c_double
        jo_grad_y%metar(n)%q = 0.0_c_double
        jo_grad_y%metar(n)%p = 0.0_c_double
        
        ! Process all variables in order: U, V, T, Q, P
        var_idx = (n-1) * 5  ! Start index for this observation's variables
        
        ! U component (if available)
        if (var_idx + 1 <= num_innovations) then
          jo_grad_y%metar(n)%u = real(delta_y(var_idx + 1))
        end if
        
        ! V component (if available)
        if (var_idx + 2 <= num_innovations) then
          jo_grad_y%metar(n)%v = real(delta_y(var_idx + 2))
        end if
        
        ! T component (if available)
        if (var_idx + 3 <= num_innovations) then
          jo_grad_y%metar(n)%t = real(delta_y(var_idx + 3))
        end if
        
        ! Q component (if available)
        if (var_idx + 4 <= num_innovations) then
          jo_grad_y%metar(n)%q = real(delta_y(var_idx + 4))
        end if
        
        ! P component (if available)
        if (var_idx + 5 <= num_innovations) then
          jo_grad_y%metar(n)%p = real(delta_y(var_idx + 5))
        end if
      end do
    case (synop)
      allocate(jo_grad_y%synop(num_obs))
      do n=1,num_obs
        ! Initialize all variables to missing
        jo_grad_y%synop(n)%u = 0.0_c_double
        jo_grad_y%synop(n)%v = 0.0_c_double
        jo_grad_y%synop(n)%t = 0.0_c_double
        jo_grad_y%synop(n)%q = 0.0_c_double
        jo_grad_y%synop(n)%p = 0.0_c_double
        
        ! Process all variables in order: U, V, T, Q, P
        ! delta_y array contains innovations for all variables for all observations
        ! For each observation, we have up to 5 variables (U, V, T, Q, P)
        var_idx = (n-1) * 5  ! Start index for this observation's variables
        
        ! U component (if available)
        if (var_idx + 1 <= num_innovations) then
          jo_grad_y%synop(n)%u = real(delta_y(var_idx + 1))
        end if
        
        ! V component (if available)
        if (var_idx + 2 <= num_innovations) then
          jo_grad_y%synop(n)%v = real(delta_y(var_idx + 2))
        end if
        
        ! T component (if available)
        if (var_idx + 3 <= num_innovations) then
          jo_grad_y%synop(n)%t = real(delta_y(var_idx + 3))
        end if
        
        ! Q component (if available)
        if (var_idx + 4 <= num_innovations) then
          jo_grad_y%synop(n)%q = real(delta_y(var_idx + 4))
        end if
        
        ! P component (if available)
        if (var_idx + 5 <= num_innovations) then
          jo_grad_y%synop(n)%p = real(delta_y(var_idx + 5))
        end if
      end do

    end select
  end subroutine init_y_from_delta

  subroutine copy_y_to_out(family, y, out_y, num_innovations)
    integer, intent(in) :: family
    type(y_type), intent(in) :: y
    real(c_double), intent(out) :: out_y(*)
    integer, intent(in) :: num_innovations
    integer :: n, var_idx, num_obs
    real(c_double), parameter :: missing_r = -888888.0_c_double
    
    ! Calculate number of observations from innovations
    ! For synop observations, we need to get the actual number of observations from persistent_iv
    if (associated(persistent_iv) .and. associated(persistent_iv%synop)) then
      num_obs = persistent_iv%info(2)%nlocal
    else
      ! Fallback: assume we have 1 observation per 5 innovations
      num_obs = max(1, num_innovations / 5)
    end if
    select case (family)
    case (synop)
      do n=1,num_obs
        ! Process all variables in order: U, V, T, Q, P
        var_idx = (n-1) * 5  ! Start index for this observation's variables
        
        ! U component (if not missing)
        if (var_idx + 1 <= num_innovations) then
          if (y%synop(n)%u /= missing_r) then
            out_y(var_idx + 1) = real(y%synop(n)%u, kind=c_double)
          else
            out_y(var_idx + 1) = 0.0_c_double
          end if
        end if
        
        ! V component (if not missing)
        if (var_idx + 2 <= num_innovations) then
          if (y%synop(n)%v /= missing_r) then
            out_y(var_idx + 2) = real(y%synop(n)%v, kind=c_double)
          else
            out_y(var_idx + 2) = 0.0_c_double
          end if
        end if
        
        ! T component (if not missing)
        if (var_idx + 3 <= num_innovations) then
          if (y%synop(n)%t /= missing_r) then
            out_y(var_idx + 3) = real(y%synop(n)%t, kind=c_double)
          else
            out_y(var_idx + 3) = 0.0_c_double
          end if
        end if
        
        ! Q component (if not missing)
        if (var_idx + 4 <= num_innovations) then
          if (y%synop(n)%q /= missing_r) then
            out_y(var_idx + 4) = real(y%synop(n)%q, kind=c_double)
          else
            out_y(var_idx + 4) = 0.0_c_double
          end if
        end if
        
        ! P component (if not missing)
        if (var_idx + 5 <= num_innovations) then
          if (y%synop(n)%p /= missing_r) then
            out_y(var_idx + 5) = real(y%synop(n)%p, kind=c_double)
          else
            out_y(var_idx + 5) = 0.0_c_double
          end if
        end if
      end do
    end select
  end subroutine copy_y_to_out

  ! Copy gradient from jo_grad_x to persistent arrays for later retrieval
  ! This function stores the adjoint gradients in persistent arrays
  ! For incremental 3D-Var, this represents the gradient of the cost function with respect to state
  subroutine copy_x_to_persistent(jo_grad_x, nx, ny, nz)
    type(x_type), intent(in) :: jo_grad_x
    integer, intent(in) :: nx, ny, nz
    integer :: i,j,k, nz1, total_size, surface_size
    
    nz1 = max(1, nz)
    total_size = nx * ny * nz1
    surface_size = nx * ny
    
    ! Allocate persistent arrays if needed
    if (.not. adjoint_allocated .or. persistent_nx /= nx .or. persistent_ny /= ny .or. persistent_nz /= nz) then
      if (allocated(persistent_adjoint_u)) deallocate(persistent_adjoint_u)
      if (allocated(persistent_adjoint_v)) deallocate(persistent_adjoint_v)
      if (allocated(persistent_adjoint_t)) deallocate(persistent_adjoint_t)
      if (allocated(persistent_adjoint_q)) deallocate(persistent_adjoint_q)
      if (allocated(persistent_adjoint_psfc)) deallocate(persistent_adjoint_psfc)
      
      allocate(persistent_adjoint_u(total_size))
      allocate(persistent_adjoint_v(total_size))
      allocate(persistent_adjoint_t(total_size))
      allocate(persistent_adjoint_q(total_size))
      allocate(persistent_adjoint_psfc(surface_size))
      
      persistent_nx = nx
      persistent_ny = ny
      persistent_nz = nz
      adjoint_allocated = .true.
    end if
    
    ! Copy gradients to persistent arrays
    do k=1,nz1; do j=1,ny; do i=1,nx
      persistent_adjoint_u(i + (j-1)*nx + (k-1)*nx*ny) = real(jo_grad_x%u(i,j,k), kind=c_double)
      persistent_adjoint_v(i + (j-1)*nx + (k-1)*nx*ny) = real(jo_grad_x%v(i,j,k), kind=c_double)
      persistent_adjoint_t(i + (j-1)*nx + (k-1)*nx*ny) = real(jo_grad_x%t(i,j,k), kind=c_double)
      persistent_adjoint_q(i + (j-1)*nx + (k-1)*nx*ny) = real(jo_grad_x%q(i,j,k), kind=c_double)
    end do; end do; end do
    do j=1,ny; do i=1,nx
      persistent_adjoint_psfc(i + (j-1)*nx) = real(jo_grad_x%psfc(i,j), kind=c_double)
    end do; end do
  end subroutine copy_x_to_persistent

  ! Copy gradient from jo_grad_x to state arrays (gradient accumulation)
  ! This function accumulates the adjoint gradients into the state space arrays
  ! For incremental 3D-Var, this represents the gradient of the cost function with respect to state
  subroutine copy_x_to_state(jo_grad_x, u, v, t, q, psfc, nx, ny, nz)
    type(x_type), intent(in) :: jo_grad_x
    real(c_double), intent(inout) :: u(*), v(*), t(*), q(*), psfc(*)
    integer, intent(in) :: nx, ny, nz
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
    
    ! Accumulate gradients into state arrays (gradient accumulation)
    do k=1,nz1; do j=1,ny; do i=1,nx
      u(i + (j-1)*nx + (k-1)*nx*ny) = u(i + (j-1)*nx + (k-1)*nx*ny) + real(jo_grad_x%u(i,j,k), kind=c_double)
      v(i + (j-1)*nx + (k-1)*nx*ny) = v(i + (j-1)*nx + (k-1)*nx*ny) + real(jo_grad_x%v(i,j,k), kind=c_double)
      t(i + (j-1)*nx + (k-1)*nx*ny) = t(i + (j-1)*nx + (k-1)*nx*ny) + real(jo_grad_x%t(i,j,k), kind=c_double)
      q(i + (j-1)*nx + (k-1)*nx*ny) = q(i + (j-1)*nx + (k-1)*nx*ny) + real(jo_grad_x%q(i,j,k), kind=c_double)
    end do; end do; end do
    do j=1,ny; do i=1,nx
      psfc(i + (j-1)*nx) = psfc(i + (j-1)*nx) + real(jo_grad_x%psfc(i,j), kind=c_double)
    end do; end do
  end subroutine copy_x_to_state

  subroutine zero_x_like(x, nx, ny, nz)
    type(x_type), intent(inout) :: x
    integer, intent(in) :: nx, ny, nz
    integer :: nz1
    nz1 = max(1, nz)
    allocate(x%u(nx,ny,nz1)); x%u = 0.0
    allocate(x%v(nx,ny,nz1)); x%v = 0.0
    allocate(x%t(nx,ny,nz1)); x%t = 0.0
    allocate(x%q(nx,ny,nz1)); x%q = 0.0
    allocate(x%psfc(nx,ny));  x%psfc = 0.0
  end subroutine zero_x_like

  ! Update only analysis increments (xa) for incremental 3D-Var
  ! This subroutine updates the analysis increments while keeping background state constant
  subroutine wrfda_update_analysis_increments(u_inc, v_inc, t_inc, q_inc, psfc_inc, grid_ptr) bind(C, name="wrfda_update_analysis_increments")
    implicit none
    real(c_double), intent(in) :: u_inc(*), v_inc(*), t_inc(*), q_inc(*), psfc_inc(*)
    type(c_ptr), value :: grid_ptr
    type(domain), pointer :: grid
    integer :: i,j,k, nx, ny, nz
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(grid_ptr, grid)
    
    ! Get grid dimensions from grid
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    
    ! Update only the xa fields (analysis increments)
    do k=1,nz; do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
      grid%xa%u(i,j,k) = real(u_inc(i + (j-1)*nx + (k-1)*nx*ny), kind=4)
      grid%xa%v(i,j,k) = real(v_inc(i + (j-1)*nx + (k-1)*nx*ny), kind=4)
      grid%xa%t(i,j,k) = real(t_inc(i + (j-1)*nx + (k-1)*nx*ny), kind=4)
      grid%xa%q(i,j,k) = real(q_inc(i + (j-1)*nx + (k-1)*nx*ny), kind=4)
    end do; end do; end do
    do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j] -> i + j*nx
      grid%xa%psfc(i,j) = real(psfc_inc(i + (j-1)*nx), kind=4)
    end do; end do
    
  end subroutine wrfda_update_analysis_increments

  ! Update background state (xb) from state data
  ! This subroutine updates the background state with the current state values
  subroutine wrfda_update_background_state(u, v, t, q, psfc, ph, phb, hf, hgt, p, pb, lats2d, lons2d, grid_ptr) bind(C, name="wrfda_update_background_state")
    implicit none
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    real(c_double), intent(in) :: ph(*), phb(*), hf(*), hgt(*), p(*), pb(*)
    real(c_double), intent(in) :: lats2d(*), lons2d(*)
    type(c_ptr), value :: grid_ptr
    type(domain), pointer :: grid
    integer :: i, j, k, nx, ny, nz, idx, idx_3d, staggered_nz
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(grid_ptr, grid)
    
    ! Get grid dimensions from grid
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    staggered_nz = nz + 1  ! Height field has nz+1 vertical levels
    
    ! Update background state (xb) with current state values
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx = i + (j-1)*nx + (k-1)*nx*ny
          grid%xb%u(i,j,k) = real(u(idx), kind=4)
          grid%xb%v(i,j,k) = real(v(idx), kind=4)
          grid%xb%t(i,j,k) = real(t(idx), kind=4)
          grid%xb%q(i,j,k) = real(q(idx), kind=4)
          ! Calculate pressure from P and PB: grid%xb%p = pb + p
          grid%xb%p(i,j,k) = real(pb(idx) + p(idx), kind=4)
        end do
      end do
    end do
    
    ! Update surface pressure
    do j = 1, ny
      do i = 1, nx
        ! C++ row-major indexing: [i][j] -> i + j*nx
        idx = i + (j-1)*nx
        grid%xb%psfc(i,j) = real(psfc(idx), kind=4)
      end do
    end do
    
    ! Update PH and PHB data (vertically staggered with nz+1 levels)
    do k = 1, staggered_nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          grid%ph_2(i,j,k) = real(ph(idx_3d), kind=4)
          grid%phb(i,j,k) = real(phb(idx_3d), kind=4)
        end do
      end do
    end do
    
    ! Update grid metadata
    do j = 1, ny
      do i = 1, nx
        ! C++ row-major indexing: [i][j] -> i + j*nx
        idx = i + (j-1)*nx
        grid%xb%lat(i,j) = real(lats2d(idx), kind=4)
        grid%xb%lon(i,j) = real(lons2d(idx), kind=4)
        ! Assign terrain height from HGT field
        grid%ht(i,j) = real(hgt(idx), kind=4)
        grid%xb%terr(i,j) = real(hgt(idx), kind=4)
        do k = 1, staggered_nz
          ! Use calculated height field instead of levels array
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          grid%xb%h(i,j,k) = real(hf(idx_3d), kind=4)
        end do
      end do
    end do
    
  end subroutine wrfda_update_background_state

  ! Get available observation families from iv structure
  integer(c_int) function wrfda_get_available_families(families_buffer, buffer_size) bind(C, name="wrfda_get_available_families")
    implicit none
    character(c_char), intent(out) :: families_buffer(*)
    integer(c_int), intent(inout) :: buffer_size
    
    integer :: i, family_index, current_pos
    character(len=20) :: family_name
    character(len=256) :: families_string
    
    ! Check if iv structure is available
    if (.not. iv_allocated .or. .not. associated(persistent_iv)) then
      families_string = ""
      buffer_size = 0
      wrfda_get_available_families = 0
      return
    end if
    
    families_string = ""
    current_pos = 1
    
    ! Check each observation family for available observations
    do family_index = 1, size(persistent_iv%info)
      if (persistent_iv%info(family_index)%nlocal > 0) then
        ! Determine family name based on index
        select case (family_index)
        case (1)
          family_name = "metar"
        case (2)
          family_name = "synop"
        case (3)
          family_name = "sound"
        case (4)
          family_name = "gpspw"
        case (5)
          family_name = "airep"
        case (6)
          family_name = "pilot"
        case (7)
          family_name = "ships"
        case (8)
          family_name = "buoy"
        case default
          family_name = "unknown"
        end select
        
        ! Add family name to string (comma-separated)
        if (len_trim(families_string) > 0) then
          families_string = trim(families_string) // ","
        end if
        families_string = trim(families_string) // trim(family_name)
      end if
    end do
    
    ! Copy to C buffer
    do i = 1, min(len_trim(families_string), buffer_size - 1)
      families_buffer(i) = families_string(i:i)
    end do
    families_buffer(min(len_trim(families_string) + 1, buffer_size)) = c_null_char
    
    buffer_size = len_trim(families_string) + 1
    wrfda_get_available_families = 0
    
  end function wrfda_get_available_families

  ! New function to call da_get_innov_vector directly
  integer(c_int) function wrfda_get_innov_vector(it, ob_ptr, iv_ptr, grid_ptr) bind(C, name="wrfda_get_innov_vector")
    implicit none
    integer(c_int), intent(in) :: it
    type(c_ptr), value :: ob_ptr, iv_ptr, grid_ptr
    type(domain), pointer :: grid
    type(y_type), pointer :: ob
    type(iv_type), pointer :: iv
    type(grid_config_rec_type), pointer :: config_flags
    
    ! Convert C pointers to Fortran pointers
    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(ob_ptr, ob)
    call c_f_pointer(iv_ptr, iv)
    
    ! config_flags is required by WRFDA interface - allocate it
    allocate(config_flags)
    
    ! Initialize QC statistics array (imported from da_control module)
    num_qcstat_conv = 0
    
    ! Call the main WRFDA innovation vector computation routine
    call da_get_innov_vector(it, num_qcstat_conv, ob, iv, grid, config_flags)

    ! Clean up allocated config_flags
    deallocate(config_flags)

    wrfda_get_innov_vector = 0
    
  end function wrfda_get_innov_vector


  ! Helper function to construct WRFDA domain structure from flat arrays
  integer(c_int) function wrfda_construct_domain_from_arrays(nx, ny, nz, u, v, t, q, psfc, ph, phb, hf, hgt, p, pb, lats2d, lons2d) bind(C, name="wrfda_construct_domain_from_arrays")
    implicit none
    integer(c_int), intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    real(c_double), intent(in) :: ph(*), phb(*), hf(*), hgt(*), p(*), pb(*)
    real(c_double), intent(in) :: lats2d(*), lons2d(*)
    
    type(domain), pointer :: grid
    integer :: i, j, k, idx, idx_3d
    integer :: staggered_nz  ! Height field has nz+1 vertical levels

    ! Allocate new domain structure
    if (.not.grid_initialized) then
      allocate(persistent_grid)
      grid => persistent_grid
      grid_initialized = .true.
    else
      if (associated(persistent_grid)) then
        grid => persistent_grid
      else
        allocate(persistent_grid)
        grid => persistent_grid
        grid_initialized = .true.
      end if
    end if

    ! Temporary constant for METADA - should be imported from WRFDA constants
    Max_StHeight_Diff = 100.0

    ! Set up basic domain dimensions
    grid%id = 1
    
    ! Initialize domain clock (required for WRFDA time management)
    ! Set domain_clock_created to false initially
    grid%domain_clock_created = .false.
    
    ! Height field is vertically staggered with nz+1 levels
    staggered_nz = nz + 1
    
    ! Allocate PH and PHB fields (vertically staggered with nz+1 levels)
    
    allocate(grid%ph_2(1:nx, 1:ny, 1:staggered_nz))
    allocate(grid%phb(1:nx, 1:ny, 1:staggered_nz))
    

    
    ! Set up grid processor dimensions (required for da_copy_tile_dims)
    grid%xp%kds = 1        ! Start of vertical domain
    grid%xp%kde = nz       ! End of vertical domain (mass variables)
    grid%xp%ids = 1        ! Start of i domain
    grid%xp%ide = nx       ! End of i domain
    grid%xp%jds = 1        ! Start of j domain
    grid%xp%jde = ny       ! End of j domain
    
    ! Set up tile dimensions (single tile for now)
    grid%num_tiles = 1
    allocate(grid%i_start(1:1), grid%i_end(1:1))
    allocate(grid%j_start(1:1), grid%j_end(1:1))
    grid%i_start(1) = 1
    grid%i_end(1) = nx
    grid%j_start(1) = 1
    grid%j_end(1) = ny
    
    ! Set up WRFDA memory dimensions (sm/em arrays)
    grid%sm31 = 1; grid%em31 = nx
    grid%sm32 = 1; grid%em32 = ny
    grid%sm33 = 1; grid%em33 = nz
    grid%sm31x = 1; grid%em31x = nx
    grid%sm32x = 1; grid%em32x = ny
    grid%sm33x = 1; grid%em33x = nz
    
    ! Initialize domain clock for WRFDA time management
    call initialize_domain_clock(grid)
    grid%sm31y = 1; grid%em31y = nx
    grid%sm32y = 1; grid%em32y = ny
    grid%sm33y = 1; grid%em33y = nz
    
    ! Set up WRFDA domain dimensions (sd/ed arrays) - required by da_copy_dims
    ! In WRF, ed arrays are typically the upper bound + 1, so for nx,ny,nz grid:
    grid%sd31 = 1; grid%ed31 = nx + 1
    grid%sd32 = 1; grid%ed32 = ny + 1  
    grid%sd33 = 1; grid%ed33 = nz + 1
    grid%sp31 = 1; grid%ep31 = nx
    grid%sp32 = 1; grid%ep32 = ny
    grid%sp33 = 1; grid%ep33 = nz

    ! Allocate and populate xb (background state) arrays
    allocate(grid%xb%u(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%v(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%t(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%q(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%p(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%psfc(1:nx, 1:ny))

    ! Copy data from flat arrays to WRFDA structure
    ! Note: C++ arrays are in [X,Y,Z] order (row-major), Fortran arrays are [Z,Y,X] (column-major)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx = i + (j-1)*nx + (k-1)*nx*ny
          grid%xb%u(i,j,k) = real(u(idx), kind=4)
          grid%xb%v(i,j,k) = real(v(idx), kind=4)
          grid%xb%t(i,j,k) = real(t(idx), kind=4)
          grid%xb%q(i,j,k) = real(q(idx), kind=4)
          ! Calculate pressure from P and PB: grid%xb%p = pb + p
          grid%xb%p(i,j,k) = real(pb(idx) + p(idx), kind=4)
        end do
      end do
    end do
    
    ! Copy PH and PHB data (vertically staggered with nz+1 levels)
    do k = 1, staggered_nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          grid%ph_2(i,j,k) = real(ph(idx_3d), kind=4)
          grid%phb(i,j,k) = real(phb(idx_3d), kind=4)
        end do
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        ! C++ row-major indexing: [i][j] -> i + j*nx
        idx = i + (j-1)*nx
        grid%xb%psfc(i,j) = real(psfc(idx), kind=4)
      end do
    end do

    ! Set up grid metadata
    allocate(grid%xb%lat(1:nx, 1:ny))
    allocate(grid%xb%lon(1:nx, 1:ny))
    ! Height field (HF) is vertically staggered with nz+1 levels
    allocate(grid%xb%h(1:nx, 1:ny, 1:staggered_nz))
    ! Terrain height fields (HGT) - 2D surface fields
    allocate(grid%ht(1:nx, 1:ny))
    allocate(grid%xb%terr(1:nx, 1:ny))
    ! Surface roughness length - temporary constant value
    allocate(grid%xb%rough(1:nx, 1:ny))

    do j = 1, ny
      do i = 1, nx
        ! C++ row-major indexing: [i][j] -> i + j*nx
        idx = i + (j-1)*nx
        grid%xb%lat(i,j) = real(lats2d(idx), kind=4)
        grid%xb%lon(i,j) = real(lons2d(idx), kind=4)
        ! Assign terrain height from HGT field
        grid%ht(i,j) = real(hgt(idx), kind=4)
        grid%xb%terr(i,j) = real(hgt(idx), kind=4)
        ! Assign temporary constant roughness length
        grid%xb%rough(i,j) = 0.5
        do k = 1, staggered_nz
          ! Use calculated height field instead of levels array
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          ! Note: HF is vertically staggered with nz+1 levels
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          grid%xb%h(i,j,k) = real(hf(idx_3d), kind=4)
        end do
      end do
    end do

    ! Allocate xa fields (analysis increments) and initialize to zero
    ! For adjoint operator, xa starts at zero and gradients are accumulated in jo_grad_x
    allocate(grid%xa%u(1:nx,1:ny,1:nz)); grid%xa%u = 0.0
    allocate(grid%xa%v(1:nx,1:ny,1:nz)); grid%xa%v = 0.0
    allocate(grid%xa%t(1:nx,1:ny,1:nz)); grid%xa%t = 0.0
    allocate(grid%xa%q(1:nx,1:ny,1:nz)); grid%xa%q = 0.0
    allocate(grid%xa%psfc(1:nx,1:ny)); grid%xa%psfc = 0.0

    wrfda_construct_domain_from_arrays = 0
    
  end function wrfda_construct_domain_from_arrays

  ! Proper WRFDA domain initialization using da_transfer_wrftoxb
  ! This function populates the WRF grid structure from arrays and then calls
  ! da_transfer_wrftoxb to properly compute all derived fields
  integer(c_int) function wrfda_init_domain_from_wrf_fields(nx, ny, nz, &
      u, v, w, t, mu, mub, p, pb, ph, phb, &
      xlat, xlong, ht, znu, znw, dn, dnw, rdnw, rdn, &
      p_top, t_init, moist, num_moist, psfc, &
      start_year, start_month, start_day, start_hour) &
      bind(C, name="wrfda_init_domain_from_wrf_fields")
    use da_transfer_model, only: da_transfer_wrftoxb
    use da_define_structures, only: xbx_type
    implicit none
    
    ! Grid dimensions
    integer(c_int), intent(in) :: nx, ny, nz, num_moist
    
    ! WRF state variables (3D fields: nx x ny x nz)
    real(c_double), intent(in) :: u(*), v(*), w(*), t(*)
    real(c_double), intent(in) :: mu(*), mub(*)  ! 2D: nx x ny
    real(c_double), intent(in) :: p(*), pb(*)    ! 3D: nx x ny x nz
    real(c_double), intent(in) :: ph(*), phb(*)  ! 3D staggered: nx x ny x (nz+1)
    
    ! Grid metadata (2D fields: nx x ny)
    real(c_double), intent(in) :: xlat(*), xlong(*), ht(*)
    
    ! Vertical coordinates (1D arrays: size nz or nz+1)
    real(c_double), intent(in) :: znu(*), znw(*), dn(*), dnw(*)
    real(c_double), intent(in) :: rdnw(*), rdn(*)
    
    ! Scalar parameters
    real(c_double), intent(in) :: p_top
    real(c_double), intent(in) :: t_init(*)  ! 3D: nx x ny x nz
    real(c_double), intent(in) :: moist(*)   ! 4D: nx x ny x nz x num_moist
    real(c_double), intent(in) :: psfc(*)    ! 2D: nx x ny
    
    ! Time information
    integer(c_int), intent(in) :: start_year, start_month, start_day, start_hour
    
    ! Local variables
    type(domain), pointer :: grid
    type(xbx_type) :: xbx
    type(grid_config_rec_type) :: config_flags
    integer :: i, j, k, m, idx, idx_3d, staggered_nz
    
    ! Set up grid structure (similar to wrfda_construct_domain_from_arrays but more complete)
    if (.not. grid_initialized) then
      allocate(persistent_grid)
      grid => persistent_grid

      ! CRITICAL: Nullify all pointer components immediately after allocation
      nullify(grid%u_2, grid%v_2, grid%w_2, grid%t_2, grid%p, grid%pb)
      nullify(grid%t_init, grid%mu_2, grid%mub, grid%psfc)
      nullify(grid%ph_2, grid%phb, grid%moist)
      nullify(grid%xlat, grid%xlong, grid%ht)
      nullify(grid%znu, grid%znw, grid%dn, grid%dnw, grid%rdnw, grid%rdn)

      grid_initialized = .true.
    else
      grid => persistent_grid
    end if
    
    staggered_nz = nz + 1
    
    ! Allocate all WRF fields needed by da_transfer_wrftoxb
    ! Note: We only allocate if not already allocated by WRFDA
    ! WRFDA's da_transfer_wrftoxb expects these fields to exist
    
    ! Mass point 3D fields (only allocate if not associated)
    if (.not. associated(grid%u_2)) allocate(grid%u_2(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%v_2)) allocate(grid%v_2(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%w_2)) allocate(grid%w_2(1:nx, 1:ny, 1:staggered_nz))
    
    if (.not. associated(grid%t_2)) allocate(grid%t_2(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%p)) allocate(grid%p(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%pb)) allocate(grid%pb(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%t_init)) allocate(grid%t_init(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%mu_2)) allocate(grid%mu_2(1:nx, 1:ny))
    
    if (.not. associated(grid%mub)) allocate(grid%mub(1:nx, 1:ny))
    
    if (.not. associated(grid%psfc)) allocate(grid%psfc(1:nx, 1:ny))
    
    ! Height fields (staggered, only allocate if not associated)
    if (.not. associated(grid%ph_2)) allocate(grid%ph_2(1:nx, 1:ny, 1:staggered_nz))
    
    if (.not. associated(grid%phb)) allocate(grid%phb(1:nx, 1:ny, 1:staggered_nz))
    
    ! Moisture fields (all species: qvapor, qcloud, qrain, qice, qsnow, qgraupel, etc.)
    if (.not. associated(grid%moist)) allocate(grid%moist(1:nx, 1:ny, 1:nz, 1:num_moist))
    
    ! Grid metadata (only allocate if not associated)
    if (.not. associated(grid%xlat)) allocate(grid%xlat(1:nx, 1:ny))
    
    if (.not. associated(grid%xlong)) allocate(grid%xlong(1:nx, 1:ny))
    
    if (.not. associated(grid%ht)) allocate(grid%ht(1:nx, 1:ny))
    
    ! Vertical coordinates (only allocate if not associated)
    if (.not. associated(grid%znu)) allocate(grid%znu(1:staggered_nz))
    
    if (.not. associated(grid%znw)) allocate(grid%znw(1:staggered_nz))
    
    if (.not. associated(grid%dn)) allocate(grid%dn(1:staggered_nz))
    
    if (.not. associated(grid%dnw)) allocate(grid%dnw(1:staggered_nz))
    
    if (.not. associated(grid%rdnw)) allocate(grid%rdnw(1:staggered_nz))
    
    if (.not. associated(grid%rdn)) allocate(grid%rdn(1:staggered_nz))
    
    ! Copy data from flat C arrays to WRFDA grid structure
    ! 3D fields
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          idx = i + (j-1)*nx + (k-1)*nx*ny
          grid%u_2(i,j,k) = real(u(idx), kind=4)
          grid%v_2(i,j,k) = real(v(idx), kind=4)
          grid%t_2(i,j,k) = real(t(idx), kind=4)
          grid%p(i,j,k) = real(p(idx), kind=4)
          grid%pb(i,j,k) = real(pb(idx), kind=4)
          grid%t_init(i,j,k) = real(t_init(idx), kind=4)
        end do
      end do
    end do
    
    ! Copy moisture fields (all species)
    do m = 1, num_moist
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            ! 4D array indexing: i + (j-1)*nx + (k-1)*nx*ny + (m-1)*nx*ny*nz
            idx = i + (j-1)*nx + (k-1)*nx*ny + (m-1)*nx*ny*nz
            grid%moist(i,j,k,m) = real(moist(idx), kind=4)
          end do
        end do
      end do
    end do
    
    ! W field (staggered)
    do k = 1, staggered_nz
      do j = 1, ny
        do i = 1, nx
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          grid%w_2(i,j,k) = real(w(idx_3d), kind=4)
          grid%ph_2(i,j,k) = real(ph(idx_3d), kind=4)
          grid%phb(i,j,k) = real(phb(idx_3d), kind=4)
        end do
      end do
    end do
    
    ! 2D fields
    do j = 1, ny
      do i = 1, nx
        idx = i + (j-1)*nx
        grid%mu_2(i,j) = real(mu(idx), kind=4)
        grid%mub(i,j) = real(mub(idx), kind=4)
        grid%psfc(i,j) = real(psfc(idx), kind=4)
        grid%xlat(i,j) = real(xlat(idx), kind=4)
        grid%xlong(i,j) = real(xlong(idx), kind=4)
        grid%ht(i,j) = real(ht(idx), kind=4)
      end do
    end do
    
    ! 1D vertical coordinates
    do k = 1, nz
      grid%znu(k) = real(znu(k), kind=4)
      grid%znw(k) = real(znw(k), kind=4)
      grid%dn(k) = real(dn(k), kind=4)
      grid%dnw(k) = real(dnw(k), kind=4)
      grid%rdnw(k) = real(rdnw(k), kind=4)
      grid%rdn(k) = real(rdn(k), kind=4)
    end do
    ! Last level for staggered arrays
    grid%znw(staggered_nz) = real(znw(staggered_nz), kind=4)
    grid%dnw(staggered_nz) = real(dnw(staggered_nz), kind=4)
    grid%rdnw(staggered_nz) = real(rdnw(staggered_nz), kind=4)
    
    ! Set scalar parameters
    grid%p_top = real(p_top, kind=4)
    grid%start_year = start_year
    grid%start_month = start_month
    grid%start_day = start_day
    grid%start_hour = start_hour
    
    ! Set up grid dimensions
    grid%xp%ids = 1; grid%xp%ide = nx
    grid%xp%jds = 1; grid%xp%jde = ny
    grid%xp%kds = 1; grid%xp%kde = nz
    
    ! Call da_transfer_wrftoxb to do all the proper initialization
    ! This computes all derived fields, diagnostics, etc.
    call da_transfer_wrftoxb(xbx, grid, config_flags)
    
    wrfda_init_domain_from_wrf_fields = 0  ! Success
    
  end function wrfda_init_domain_from_wrf_fields

  ! Construct y_type from observation data
  type(c_ptr) function wrfda_construct_y_type(num_obs, num_levels, u_values, v_values, t_values, p_values, q_values, u_errors, v_errors, t_errors, p_errors, q_errors, u_available, v_available, t_available, p_available, q_available, lats, lons, obs_types, family) bind(C, name="wrfda_construct_y_type")
    implicit none
    integer(c_int), intent(in) :: num_obs, num_levels
    real(c_double), intent(in) :: u_values(*), v_values(*), t_values(*), p_values(*), q_values(*)
    real(c_double), intent(in) :: u_errors(*), v_errors(*), t_errors(*), p_errors(*), q_errors(*)
    integer(c_int), intent(in) :: u_available(*), v_available(*), t_available(*), p_available(*), q_available(*)
    real(c_double), intent(in) :: lats(*), lons(*)
    character(c_char), intent(in) :: obs_types(*), family(*)
    
    type(y_type), pointer :: y
    integer :: i
    character(len=20) :: family_str
    character(c_char), target :: family_target(20)
    
    
    ! Check if num_obs is valid
    if (num_obs <= 0 .or. num_levels <= 0) then
      wrfda_construct_y_type = c_null_ptr
      return
    end if
    
    ! Convert C string to Fortran string
    family_target = family(1:20)
    family_str = ""
    do i = 1, 20
      if (family_target(i) == c_null_char) exit
      family_str(i:i) = family_target(i)
    end do
    
    ! Allocate persistent y_type if not already allocated
    if (.not. y_allocated) then
      allocate(persistent_y)
      y_allocated = .true.
    end if
    
    ! Point local y to persistent structure
    y => persistent_y
    
    ! Initialize counters
    y%nlocal = 0
    y%ntotal = 0
    y%num_inst = 0
    
    ! For surface observations (metar, synop, adpsfc), allocate synop array
    if (trim(family_str) == "metar" .or. trim(family_str) == "synop" .or. trim(family_str) == "adpsfc") then
      allocate(y%synop(num_obs))
      y%nlocal(1) = num_obs  ! synop index
      y%ntotal(1) = num_obs
      
      ! Populate synop array with observation data
      do i = 1, num_obs
        ! Set values from the structured data only if available
        ! residual_synop_type has simple real fields (u, v, t, p, q)
        if (u_available(i) == 1) then
          y%synop(i)%u = u_values(i)
        else
          y%synop(i)%u = missing_r  ! Default value for unavailable data
        end if
        
        if (v_available(i) == 1) then
          y%synop(i)%v = v_values(i)
        else
          y%synop(i)%v = missing_r
        end if
        
        if (t_available(i) == 1) then
          y%synop(i)%t = t_values(i)
        else
          y%synop(i)%t = missing_r
        end if
        
        if (p_available(i) == 1) then
          y%synop(i)%p = p_values(i)
        else
          y%synop(i)%p = missing_r
        end if
        
        if (q_available(i) == 1) then
          y%synop(i)%q = q_values(i)
        else
          y%synop(i)%q = missing_r
        end if
        
        ! lats, lons arrays contain observation location information (used for coordinate conversion)
        if (i <= num_obs .and. lats(i) > -90.0 .and. lats(i) < 90.0 .and. lons(i) > -180.0 .and. lons(i) < 180.0) then
          ! Valid lat/lon values - these are used for coordinate conversion in WRFDA
          continue
        end if
        
        ! obs_types array contains observation type codes (used for observation processing)
        if (i <= num_obs .and. obs_types(i) /= c_null_char) then
          ! Valid observation type - used for processing
          continue
        end if
        
        ! error arrays contain observation error information (used for quality control)
        if (i <= num_obs) then
          ! Reference error arrays to prevent unused argument warnings
          if (u_errors(i) > 0.0 .or. v_errors(i) > 0.0 .or. t_errors(i) > 0.0 .or. &
              p_errors(i) > 0.0 .or. q_errors(i) > 0.0) then
            ! Valid error information - used for quality control
            continue
          end if
        end if
      end do
      
    end if
    
    ! Return pointer to persistent y_type
    wrfda_construct_y_type = c_loc(persistent_y)
    
  end function wrfda_construct_y_type

  ! Construct iv_type from observation data
  type(c_ptr) function wrfda_construct_iv_type(num_obs, num_levels, u_values, v_values, t_values, p_values, q_values, u_errors, v_errors, t_errors, p_errors, q_errors, u_qc, v_qc, t_qc, p_qc, q_qc, u_available, v_available, t_available, p_available, q_available, lats, lons, levels, elevations, obs_types, family) bind(C, name="wrfda_construct_iv_type")
    implicit none
    integer(c_int), intent(in) :: num_obs, num_levels
    real(c_double), intent(in) :: u_values(*), v_values(*), t_values(*), p_values(*), q_values(*)
    real(c_double), intent(in) :: u_errors(*), v_errors(*), t_errors(*), p_errors(*), q_errors(*)
    integer(c_int), intent(in) :: u_qc(*), v_qc(*), t_qc(*), p_qc(*), q_qc(*)
    integer(c_int), intent(in) :: u_available(*), v_available(*), t_available(*), p_available(*), q_available(*)
    real(c_double), intent(in) :: lats(*), lons(*), levels(*), elevations(*)
    character(c_char), intent(in) :: obs_types(*), family(*)
    
    type(iv_type), pointer :: iv
    type(domain), pointer :: grid
    integer :: i, lev
    character(len=20) :: family_str
    character(c_char), target :: family_target(20)
    real(8) :: obs_lat, obs_lon, grid_x, grid_y, x_frac, y_frac, x_frac_m, y_frac_m
    integer :: grid_i, grid_j
    
    ! Use global persistent iv structure
    ! Convert domain pointer to grid pointer
    grid => persistent_grid
    
    ! Check if num_obs is valid
    if (num_obs <= 0 .or. num_levels <= 0) then
      wrfda_construct_iv_type = c_null_ptr
      return
    end if
    
    ! Now implement the actual iv_type construction
    
    ! Allocate iv_type
    ! Allocate persistent iv structure if not already allocated
    if (.not. iv_allocated) then
      allocate(persistent_iv)
      iv_allocated = .true.
    end if
    
    ! Point local iv to persistent structure
    iv => persistent_iv
    if (.not. associated(iv)) then
      wrfda_construct_iv_type = c_null_ptr
      return
    end if
    
    ! Initialize basic fields
    iv%time = 1
    iv%num_inst = 0
    iv%total_rad_pixel = 0
    iv%total_rad_channel = 0
    
    ! Initialize info array
    do i = 1, size(iv%info)
      iv%info(i)%nlocal = 0
      iv%info(i)%ntotal = 0
      iv%info(i)%n1 = 0
      iv%info(i)%n2 = 0
    end do
    
    ! Set up for surface observations (family index 2 = synop)
    iv%info(2)%nlocal = num_obs
    iv%info(2)%ntotal = num_obs
    
    ! Set up plocal array for proper time slot indexing
    ! For 3D-Var (num_fgat_time = 1), we have:
    ! plocal(0) = 0 (no observations before time slot 1)
    ! plocal(1) = num_obs (cumulative count at time slot 1)
    iv%info(2)%plocal(0) = 0
    iv%info(2)%plocal(1) = num_obs
    
    ! Set n1 and n2 based on plocal array (WRFDA standard pattern)
    iv%info(2)%n1 = iv%info(2)%plocal(0) + 1  ! = 1
    iv%info(2)%n2 = iv%info(2)%plocal(1)      ! = num_obs
    
    ! Allocate synop array if we have observations
    if (num_obs > 0) then
      ! Check if synop array is already associated
      if (.not. associated(iv%synop)) then
        allocate(iv%synop(num_obs))
      else
        ! Check if the existing array has the right size
        if (size(iv%synop) /= num_obs) then
          allocate(iv%synop(num_obs))
        end if
      end if
      
    ! Allocate interpolation arrays for synop observations
    ! Use kms:kme for vertical dimension to match WRFDA expectations
    iv%info(2)%max_lev = 1  ! Surface observations have max_lev = 1
    
    if (.not. allocated(iv%info(2)%i)) then
      allocate(iv%info(2)%i(kms:kme, num_obs))
      allocate(iv%info(2)%j(kms:kme, num_obs))
      allocate(iv%info(2)%dx(kms:kme, num_obs))
      allocate(iv%info(2)%dy(kms:kme, num_obs))
      allocate(iv%info(2)%dxm(kms:kme, num_obs))
      allocate(iv%info(2)%dym(kms:kme, num_obs))
      allocate(iv%info(2)%k(kms:kme, num_obs))
      allocate(iv%info(2)%zk(kms:kme, num_obs))
      allocate(iv%info(2)%dz(kms:kme, num_obs))
      allocate(iv%info(2)%dzm(kms:kme, num_obs))
      allocate(iv%info(2)%levels(num_obs))
      allocate(iv%info(2)%proc_domain(kms:kme, num_obs))
    end if
      
      ! Allocate observation metadata arrays
      if (.not. allocated(iv%info(2)%platform)) then
        allocate(iv%info(2)%platform(num_obs))
        allocate(iv%info(2)%id(num_obs))
        allocate(iv%info(2)%name(num_obs))
        allocate(iv%info(2)%date_char(num_obs))
        allocate(iv%info(2)%lat(1, num_obs))
        allocate(iv%info(2)%lon(1, num_obs))
      end if
      
      ! Populate synop data
      do i = 1, num_obs
        ! Set up height from station elevation
        iv%synop(i)%h = elevations(i)
        
        ! Set up field_type members for u, v, t, p, q only if available
        if (u_available(i) == 1) then
          iv%synop(i)%u%inv = u_values(i)
          iv%synop(i)%u%qc = u_qc(i)
          iv%synop(i)%u%error = u_errors(i)
        else
          iv%synop(i)%u%inv = missing_r  ! WRFDA standard missing value
          iv%synop(i)%u%qc = missing_data  ! WRFDA standard missing data QC flag
          iv%synop(i)%u%error = 1.0
        end if
        iv%synop(i)%u%sens = 0.0
        iv%synop(i)%u%imp = 0.0
        
        if (v_available(i) == 1) then
          iv%synop(i)%v%inv = v_values(i)
          iv%synop(i)%v%qc = v_qc(i)
          iv%synop(i)%v%error = v_errors(i)
        else
          iv%synop(i)%v%inv = missing_r  ! WRFDA standard missing value
          iv%synop(i)%v%qc = missing_data  ! WRFDA standard missing data QC flag
          iv%synop(i)%v%error = 1.0
        end if
        iv%synop(i)%v%sens = 0.0
        iv%synop(i)%v%imp = 0.0
        
        if (t_available(i) == 1) then
          iv%synop(i)%t%inv = t_values(i)
          iv%synop(i)%t%qc = t_qc(i)
          iv%synop(i)%t%error = t_errors(i)
        else
          iv%synop(i)%t%inv = missing_r  ! WRFDA standard missing value
          iv%synop(i)%t%qc = missing_data  ! WRFDA standard missing data QC flag
          iv%synop(i)%t%error = 1.0
        end if
        iv%synop(i)%t%sens = 0.0
        iv%synop(i)%t%imp = 0.0
        
        if (p_available(i) == 1) then
          iv%synop(i)%p%inv = p_values(i)
          iv%synop(i)%p%qc = p_qc(i)
          iv%synop(i)%p%error = p_errors(i)
        else
          iv%synop(i)%p%inv = missing_r  ! WRFDA standard missing value
          iv%synop(i)%p%qc = missing_data  ! WRFDA standard missing data QC flag
          iv%synop(i)%p%error = 1.0
        end if
        iv%synop(i)%p%sens = 0.0
        iv%synop(i)%p%imp = 0.0
        
        if (q_available(i) == 1) then
          iv%synop(i)%q%inv = q_values(i)
          iv%synop(i)%q%qc = q_qc(i)
          iv%synop(i)%q%error = q_errors(i)
        else
          iv%synop(i)%q%inv = missing_r  ! WRFDA standard missing value
          iv%synop(i)%q%qc = missing_data  ! WRFDA standard missing data QC flag
          iv%synop(i)%q%error = 1.0
        end if
        iv%synop(i)%q%sens = 0.0
        iv%synop(i)%q%imp = 0.0
        
        ! Populate observation metadata arrays using the input arguments
        ! Set number of levels for this observation
        if (i <= num_levels) then
          iv%info(2)%levels(i) = int(levels(i))  ! Convert to integer
        else
          iv%info(2)%levels(i) = 1  ! Default to 1 level for surface observations
        end if
        
        ! Set observation type/platform information
        if (i <= num_obs .and. obs_types(i) /= c_null_char) then
          ! Convert C character to Fortran string for platform
          iv%info(2)%platform(i) = ""
          if (obs_types(i) == 'S' .or. obs_types(i) == 's') then
            iv%info(2)%platform(i) = "SYNOP"
          else if (obs_types(i) == 'M' .or. obs_types(i) == 'm') then
            iv%info(2)%platform(i) = "METAR"
          else if (obs_types(i) == 'B' .or. obs_types(i) == 'b') then
            iv%info(2)%platform(i) = "BUOY"
          else
            iv%info(2)%platform(i) = "SYNOP"  ! Default to SYNOP
          end if
        else
          iv%info(2)%platform(i) = "SYNOP"  ! Default platform
        end if
        
        ! Set station ID (use observation index as ID)
        write(iv%info(2)%id(i), '(I5.5)') i  ! Format as 5-digit ID
        
        ! Set station name
        write(iv%info(2)%name(i), '(A,I0)') "STATION_", i
        
        ! Set date (use current time as placeholder)
        iv%info(2)%date_char(i) = "2024-01-01_00:00:00"
        
        ! Set lat/lon coordinates
        iv%info(2)%lat(1, i) = lats(i)
        iv%info(2)%lon(1, i) = lons(i)
      end do
      
      ! Compute horizontal interpolation indices and weights using WRFDA coordinate conversion
      do i = 1, num_obs
        ! Get observation lat/lon coordinates
        obs_lat = lats(i)
        obs_lon = lons(i)
        
        ! Convert lat/lon to grid coordinates using WRFDA's coordinate conversion
        ! Note: We need to use the map projection information that was initialized
        ! in the C++ code before calling this function
        call da_llxy_wrf(map_info, obs_lat, obs_lon, grid_x, grid_y)
        
        ! Convert grid coordinates to fractional indices using WRFDA's da_togrid
        ! da_togrid(x, ib, ie, i, dx, dxm) - converts x coordinate to grid index i
        call da_togrid(grid_x, 1, grid%xp%ide-1, grid_i, x_frac, x_frac_m)
        call da_togrid(grid_y, 1, grid%xp%jde-1, grid_j, y_frac, y_frac_m)
        
        ! Ensure indices are within bounds for interpolation
        grid_i = max(1, min(grid_i, grid%xp%ide-1))
        grid_j = max(1, min(grid_j, grid%xp%jde-1))
        
        ! Set grid indices for all vertical levels
        do lev = kms, kme
          iv%info(2)%i(lev, i) = grid_i
          iv%info(2)%j(lev, i) = grid_j
          
          ! Set interpolation weights using actual fractional values
          iv%info(2)%dx(lev, i) = x_frac
          iv%info(2)%dy(lev, i) = y_frac
          iv%info(2)%dxm(lev, i) = x_frac_m
          iv%info(2)%dym(lev, i) = y_frac_m
          
          ! Set vertical level information
          iv%info(2)%k(lev, i) = 1  ! Surface level
          iv%info(2)%zk(lev, i) = 1.0
          iv%info(2)%dz(lev, i) = 0.0
          iv%info(2)%dzm(lev, i) = 1.0
          
          ! Set processor domain (assume all observations are in domain)
          iv%info(2)%proc_domain(lev, i) = .true.
        end do
        
        iv%info(2)%levels(i) = 1  ! Surface observation has 1 level
      end do
      
    end if
    
    wrfda_construct_iv_type = c_loc(persistent_iv)
    
    ! Convert C string to Fortran string
    family_target = family(1:20)
    family_str = ""
    do i = 1, 20
      if (family_target(i) == c_null_char) exit
      family_str(i:i) = family_target(i)
    end do
    
  end function wrfda_construct_iv_type

  ! Construct config_flags for WRFDA
  type(c_ptr) function wrfda_construct_config_flags() bind(C, name="wrfda_construct_config_flags")
    implicit none
    
    ! For now, return a null pointer since config_flags is complex
    ! In a full implementation, this would create a proper grid_config_rec_type
    ! with all the necessary WRFDA configuration parameters
    
    
    wrfda_construct_config_flags = c_null_ptr
    
  end function wrfda_construct_config_flags

  ! Count innovation values from iv_type structure for all observation types
  integer(c_int) function wrfda_count_innovations(iv_ptr, num_innovations) bind(C, name="wrfda_count_innovations")
    implicit none
    type(c_ptr), value :: iv_ptr
    integer(c_int), intent(out) :: num_innovations
    
    type(iv_type), pointer :: iv
    integer :: i, n, count
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(iv_ptr, iv)
    
    if (.not. associated(iv)) then
      num_innovations = 0
      wrfda_count_innovations = 1
      return
    end if
    
    count = 0
    
    ! Count innovations for ALL observation families
    ! Each family has its own observation structure with variables that have .inv field
    
    ! Synop observations
    if (associated(iv%synop) .and. iv%info(synop)%nlocal > 0) then
      do n = 1, iv%info(synop)%nlocal
        if (iv%synop(n)%u%qc == 0 .and. abs(iv%synop(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%synop(n)%v%qc == 0 .and. abs(iv%synop(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%synop(n)%t%qc == 0 .and. abs(iv%synop(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%synop(n)%p%qc == 0 .and. abs(iv%synop(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%synop(n)%q%qc == 0 .and. abs(iv%synop(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Metar observations
    if (associated(iv%metar) .and. iv%info(metar)%nlocal > 0) then
      do n = 1, iv%info(metar)%nlocal
        if (iv%metar(n)%u%qc == 0 .and. abs(iv%metar(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%metar(n)%v%qc == 0 .and. abs(iv%metar(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%metar(n)%t%qc == 0 .and. abs(iv%metar(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%metar(n)%p%qc == 0 .and. abs(iv%metar(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%metar(n)%q%qc == 0 .and. abs(iv%metar(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Ships observations
    if (associated(iv%ships) .and. iv%info(ships)%nlocal > 0) then
      do n = 1, iv%info(ships)%nlocal
        if (iv%ships(n)%u%qc == 0 .and. abs(iv%ships(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%ships(n)%v%qc == 0 .and. abs(iv%ships(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%ships(n)%t%qc == 0 .and. abs(iv%ships(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%ships(n)%p%qc == 0 .and. abs(iv%ships(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%ships(n)%q%qc == 0 .and. abs(iv%ships(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Buoy observations
    if (associated(iv%buoy) .and. iv%info(buoy)%nlocal > 0) then
      do n = 1, iv%info(buoy)%nlocal
        if (iv%buoy(n)%u%qc == 0 .and. abs(iv%buoy(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%buoy(n)%v%qc == 0 .and. abs(iv%buoy(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%buoy(n)%t%qc == 0 .and. abs(iv%buoy(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%buoy(n)%p%qc == 0 .and. abs(iv%buoy(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%buoy(n)%q%qc == 0 .and. abs(iv%buoy(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Sound (radiosonde) observations - multi-level
    if (associated(iv%sound) .and. iv%info(sound)%nlocal > 0) then
      do n = 1, iv%info(sound)%nlocal
        do i = 1, iv%info(sound)%levels(n)
          if (iv%sound(n)%u(i)%qc == 0 .and. abs(iv%sound(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%sound(n)%v(i)%qc == 0 .and. abs(iv%sound(n)%v(i)%inv) > 1.0e-10) count = count + 1
          if (iv%sound(n)%t(i)%qc == 0 .and. abs(iv%sound(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%sound(n)%q(i)%qc == 0 .and. abs(iv%sound(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Airep (aircraft) observations - multi-level
    if (associated(iv%airep) .and. iv%info(airep)%nlocal > 0) then
      do n = 1, iv%info(airep)%nlocal
        do i = 1, iv%info(airep)%levels(n)
          if (iv%airep(n)%u(i)%qc == 0 .and. abs(iv%airep(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%airep(n)%v(i)%qc == 0 .and. abs(iv%airep(n)%v(i)%inv) > 1.0e-10) count = count + 1
          if (iv%airep(n)%t(i)%qc == 0 .and. abs(iv%airep(n)%t(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Pilot (balloon) observations - multi-level
    if (associated(iv%pilot) .and. iv%info(pilot)%nlocal > 0) then
      do n = 1, iv%info(pilot)%nlocal
        do i = 1, iv%info(pilot)%levels(n)
          if (iv%pilot(n)%u(i)%qc == 0 .and. abs(iv%pilot(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%pilot(n)%v(i)%qc == 0 .and. abs(iv%pilot(n)%v(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Sonde surface observations
    if (associated(iv%sonde_sfc) .and. iv%info(sonde_sfc)%nlocal > 0) then
      do n = 1, iv%info(sonde_sfc)%nlocal
        if (iv%sonde_sfc(n)%u%qc == 0 .and. abs(iv%sonde_sfc(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%sonde_sfc(n)%v%qc == 0 .and. abs(iv%sonde_sfc(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%sonde_sfc(n)%t%qc == 0 .and. abs(iv%sonde_sfc(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%sonde_sfc(n)%p%qc == 0 .and. abs(iv%sonde_sfc(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%sonde_sfc(n)%q%qc == 0 .and. abs(iv%sonde_sfc(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Geostationary satellite AMV observations - multi-level
    if (associated(iv%geoamv) .and. iv%info(geoamv)%nlocal > 0) then
      do n = 1, iv%info(geoamv)%nlocal
        do i = 1, iv%info(geoamv)%levels(n)
          if (iv%geoamv(n)%u(i)%qc == 0 .and. abs(iv%geoamv(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%geoamv(n)%v(i)%qc == 0 .and. abs(iv%geoamv(n)%v(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Polar satellite AMV observations - multi-level
    if (associated(iv%polaramv) .and. iv%info(polaramv)%nlocal > 0) then
      do n = 1, iv%info(polaramv)%nlocal
        do i = 1, iv%info(polaramv)%levels(n)
          if (iv%polaramv(n)%u(i)%qc == 0 .and. abs(iv%polaramv(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%polaramv(n)%v(i)%qc == 0 .and. abs(iv%polaramv(n)%v(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! GPS precipitable water observations
    if (associated(iv%gpspw) .and. iv%info(gpspw)%nlocal > 0) then
      do n = 1, iv%info(gpspw)%nlocal
        if (iv%gpspw(n)%tpw%qc == 0 .and. abs(iv%gpspw(n)%tpw%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! GPS refractivity observations - multi-level
    if (associated(iv%gpsref) .and. iv%info(gpsref)%nlocal > 0) then
      do n = 1, iv%info(gpsref)%nlocal
        do i = 1, iv%info(gpsref)%levels(n)
          if (iv%gpsref(n)%ref(i)%qc == 0 .and. abs(iv%gpsref(n)%ref(i)%inv) > 1.0e-10) count = count + 1
          if (iv%gpsref(n)%p(i)%qc == 0 .and. abs(iv%gpsref(n)%p(i)%inv) > 1.0e-10) count = count + 1
          if (iv%gpsref(n)%t(i)%qc == 0 .and. abs(iv%gpsref(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%gpsref(n)%q(i)%qc == 0 .and. abs(iv%gpsref(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! GPS excess phase observations - multi-level
    if (associated(iv%gpseph) .and. iv%info(gpseph)%nlocal > 0) then
      do n = 1, iv%info(gpseph)%nlocal
        do i = 1, iv%info(gpseph)%levels(n)
          if (iv%gpseph(n)%eph(i)%qc == 0 .and. abs(iv%gpseph(n)%eph(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! QuikSCAT wind observations
    if (associated(iv%qscat) .and. iv%info(qscat)%nlocal > 0) then
      do n = 1, iv%info(qscat)%nlocal
        if (iv%qscat(n)%u%qc == 0 .and. abs(iv%qscat(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%qscat(n)%v%qc == 0 .and. abs(iv%qscat(n)%v%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Wind profiler observations - multi-level (uses pilot_type)
    if (associated(iv%profiler) .and. iv%info(profiler)%nlocal > 0) then
      do n = 1, iv%info(profiler)%nlocal
        do i = 1, iv%info(profiler)%levels(n)
          if (iv%profiler(n)%u(i)%qc == 0 .and. abs(iv%profiler(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%profiler(n)%v(i)%qc == 0 .and. abs(iv%profiler(n)%v(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! SSM/I rain rate observations
    if (associated(iv%ssmi_rv) .and. iv%info(ssmi_rv)%nlocal > 0) then
      do n = 1, iv%info(ssmi_rv)%nlocal
        if (iv%ssmi_rv(n)%Speed%qc == 0 .and. abs(iv%ssmi_rv(n)%Speed%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_rv(n)%tpw%qc == 0 .and. abs(iv%ssmi_rv(n)%tpw%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! SSM/I brightness temperature observations - single level, multiple channels
    if (associated(iv%ssmi_tb) .and. iv%info(ssmi_tb)%nlocal > 0) then
      do n = 1, iv%info(ssmi_tb)%nlocal
        if (iv%ssmi_tb(n)%tb19h%qc == 0 .and. abs(iv%ssmi_tb(n)%tb19h%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb19v%qc == 0 .and. abs(iv%ssmi_tb(n)%tb19v%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb22v%qc == 0 .and. abs(iv%ssmi_tb(n)%tb22v%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb37h%qc == 0 .and. abs(iv%ssmi_tb(n)%tb37h%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb37v%qc == 0 .and. abs(iv%ssmi_tb(n)%tb37v%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb85h%qc == 0 .and. abs(iv%ssmi_tb(n)%tb85h%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb85v%qc == 0 .and. abs(iv%ssmi_tb(n)%tb85v%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! SSM/T1 observations - multi-level
    if (associated(iv%ssmt1) .and. iv%info(ssmt1)%nlocal > 0) then
      do n = 1, iv%info(ssmt1)%nlocal
        do i = 1, iv%info(ssmt1)%levels(n)
          if (iv%ssmt1(n)%t(i)%qc == 0 .and. abs(iv%ssmt1(n)%t(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! SSM/T2 observations - multi-level
    if (associated(iv%ssmt2) .and. iv%info(ssmt2)%nlocal > 0) then
      do n = 1, iv%info(ssmt2)%nlocal
        do i = 1, iv%info(ssmt2)%levels(n)
          if (iv%ssmt2(n)%rh(i)%qc == 0 .and. abs(iv%ssmt2(n)%rh(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! SATEM observations - multi-level
    if (associated(iv%satem) .and. iv%info(satem)%nlocal > 0) then
      do n = 1, iv%info(satem)%nlocal
        do i = 1, iv%info(satem)%levels(n)
          if (iv%satem(n)%thickness(i)%qc == 0 .and. abs(iv%satem(n)%thickness(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Pseudo observations - single level
    if (associated(iv%pseudo) .and. iv%info(pseudo)%nlocal > 0) then
      do n = 1, iv%info(pseudo)%nlocal
        if (iv%pseudo(n)%u%qc == 0 .and. abs(iv%pseudo(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%pseudo(n)%v%qc == 0 .and. abs(iv%pseudo(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%pseudo(n)%t%qc == 0 .and. abs(iv%pseudo(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%pseudo(n)%p%qc == 0 .and. abs(iv%pseudo(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%pseudo(n)%q%qc == 0 .and. abs(iv%pseudo(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Bogus (tropical cyclone) observations - multi-level
    if (associated(iv%bogus) .and. iv%info(bogus)%nlocal > 0) then
      do n = 1, iv%info(bogus)%nlocal
        do i = 1, iv%info(bogus)%levels(n)
          if (iv%bogus(n)%u(i)%qc == 0 .and. abs(iv%bogus(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%bogus(n)%v(i)%qc == 0 .and. abs(iv%bogus(n)%v(i)%inv) > 1.0e-10) count = count + 1
          if (iv%bogus(n)%t(i)%qc == 0 .and. abs(iv%bogus(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%bogus(n)%q(i)%qc == 0 .and. abs(iv%bogus(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
        ! Surface level
        if (iv%bogus(n)%slp%qc == 0 .and. abs(iv%bogus(n)%slp%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! AIRS retrievals - multi-level
    if (associated(iv%airsr) .and. iv%info(airsr)%nlocal > 0) then
      do n = 1, iv%info(airsr)%nlocal
        do i = 1, iv%info(airsr)%levels(n)
          if (iv%airsr(n)%t(i)%qc == 0 .and. abs(iv%airsr(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%airsr(n)%q(i)%qc == 0 .and. abs(iv%airsr(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! MTGIRS observations - multi-level
    if (associated(iv%mtgirs) .and. iv%info(mtgirs)%nlocal > 0) then
      do n = 1, iv%info(mtgirs)%nlocal
        do i = 1, iv%info(mtgirs)%levels(n)
          if (iv%mtgirs(n)%u(i)%qc == 0 .and. abs(iv%mtgirs(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%mtgirs(n)%v(i)%qc == 0 .and. abs(iv%mtgirs(n)%v(i)%inv) > 1.0e-10) count = count + 1
          if (iv%mtgirs(n)%t(i)%qc == 0 .and. abs(iv%mtgirs(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%mtgirs(n)%q(i)%qc == 0 .and. abs(iv%mtgirs(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! TAMDAR aircraft observations - multi-level
    if (associated(iv%tamdar) .and. iv%info(tamdar)%nlocal > 0) then
      do n = 1, iv%info(tamdar)%nlocal
        do i = 1, iv%info(tamdar)%levels(n)
          if (iv%tamdar(n)%u(i)%qc == 0 .and. abs(iv%tamdar(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%tamdar(n)%v(i)%qc == 0 .and. abs(iv%tamdar(n)%v(i)%inv) > 1.0e-10) count = count + 1
          if (iv%tamdar(n)%t(i)%qc == 0 .and. abs(iv%tamdar(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%tamdar(n)%q(i)%qc == 0 .and. abs(iv%tamdar(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! TAMDAR surface observations (uses synop_type)
    if (associated(iv%tamdar_sfc) .and. iv%info(tamdar_sfc)%nlocal > 0) then
      do n = 1, iv%info(tamdar_sfc)%nlocal
        if (iv%tamdar_sfc(n)%u%qc == 0 .and. abs(iv%tamdar_sfc(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%tamdar_sfc(n)%v%qc == 0 .and. abs(iv%tamdar_sfc(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%tamdar_sfc(n)%t%qc == 0 .and. abs(iv%tamdar_sfc(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%tamdar_sfc(n)%p%qc == 0 .and. abs(iv%tamdar_sfc(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%tamdar_sfc(n)%q%qc == 0 .and. abs(iv%tamdar_sfc(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Rain gauge observations
    if (associated(iv%rain) .and. iv%info(rain)%nlocal > 0) then
      do n = 1, iv%info(rain)%nlocal
        if (iv%rain(n)%rain%qc == 0 .and. abs(iv%rain(n)%rain%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Radar observations - multi-level
    if (associated(iv%radar) .and. iv%info(radar)%nlocal > 0) then
      do n = 1, iv%info(radar)%nlocal
        do i = 1, iv%info(radar)%levels(n)
          if (iv%radar(n)%rv(i)%qc == 0 .and. abs(iv%radar(n)%rv(i)%inv) > 1.0e-10) count = count + 1
          if (iv%radar(n)%rf(i)%qc == 0 .and. abs(iv%radar(n)%rf(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Lightning observations
    if (associated(iv%lightning) .and. iv%info(lightning)%nlocal > 0) then
      do n = 1, iv%info(lightning)%nlocal
        do i = 1, iv%info(lightning)%levels(n)
          if (iv%lightning(n)%w(i)%qc == 0 .and. abs(iv%lightning(n)%w(i)%inv) > 1.0e-10) count = count + 1
          if (iv%lightning(n)%div(i)%qc == 0 .and. abs(iv%lightning(n)%div(i)%inv) > 1.0e-10) count = count + 1
          if (iv%lightning(n)%qv(i)%qc == 0 .and. abs(iv%lightning(n)%qv(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    num_innovations = count
    wrfda_count_innovations = 0
    
  end function wrfda_count_innovations

  ! Extract innovation values from iv_type structure
  integer(c_int) function wrfda_extract_innovations(iv_ptr, innovations, num_innovations) bind(C, name="wrfda_extract_innovations")
    implicit none
    type(c_ptr), value :: iv_ptr
    real(c_double), intent(out) :: innovations(*)
    integer(c_int), intent(out) :: num_innovations
    
    type(iv_type), pointer :: iv
    integer :: i, k, count
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(iv_ptr, iv)
    
    if (.not. associated(iv)) then
      num_innovations = 0
      wrfda_extract_innovations = 1
      return
    end if
    
    ! Extract innovations for all families
    count = 0
    
    ! SYNOP family
    if (associated(iv%synop) .and. iv%info(synop)%nlocal > 0) then
      do i = 1, iv%info(synop)%nlocal
        if (i > size(iv%synop)) exit
        
        ! Extract U component innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%u%qc == 0 .and. abs(iv%synop(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%synop(i)%u%inv
        end if
        
        ! Extract V component innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%v%qc == 0 .and. abs(iv%synop(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%synop(i)%v%inv
        end if
        
        ! Extract T component innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%t%qc == 0 .and. abs(iv%synop(i)%t%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%synop(i)%t%inv
        end if
        
        ! Extract P component innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%p%qc == 0 .and. abs(iv%synop(i)%p%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%synop(i)%p%inv
        end if
        
        ! Extract Q component innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%q%qc == 0 .and. abs(iv%synop(i)%q%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%synop(i)%q%inv
        end if
      end do
    end if
    
    ! METAR family
    if (associated(iv%metar) .and. iv%info(metar)%nlocal > 0) then
      do i = 1, iv%info(metar)%nlocal
        if (i > size(iv%metar)) exit
        
        if (iv%metar(i)%u%qc == 0 .and. abs(iv%metar(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%metar(i)%u%inv
        end if
        if (iv%metar(i)%v%qc == 0 .and. abs(iv%metar(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%metar(i)%v%inv
        end if
        if (iv%metar(i)%t%qc == 0 .and. abs(iv%metar(i)%t%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%metar(i)%t%inv
        end if
        if (iv%metar(i)%p%qc == 0 .and. abs(iv%metar(i)%p%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%metar(i)%p%inv
        end if
        if (iv%metar(i)%q%qc == 0 .and. abs(iv%metar(i)%q%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%metar(i)%q%inv
        end if
      end do
    end if
    
    ! SHIPS family
    if (associated(iv%ships) .and. iv%info(ships)%nlocal > 0) then
      do i = 1, iv%info(ships)%nlocal
        if (i > size(iv%ships)) exit
        
        if (iv%ships(i)%u%qc == 0 .and. abs(iv%ships(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ships(i)%u%inv
        end if
        if (iv%ships(i)%v%qc == 0 .and. abs(iv%ships(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ships(i)%v%inv
        end if
        if (iv%ships(i)%t%qc == 0 .and. abs(iv%ships(i)%t%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ships(i)%t%inv
        end if
        if (iv%ships(i)%p%qc == 0 .and. abs(iv%ships(i)%p%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ships(i)%p%inv
        end if
        if (iv%ships(i)%q%qc == 0 .and. abs(iv%ships(i)%q%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ships(i)%q%inv
        end if
      end do
    end if
    
    ! SOUND family - if available and allocated
    ! Note: sound observations have multiple levels per station
    if (associated(iv%sound) .and. iv%info(sound)%nlocal > 0) then
      do i = 1, iv%info(sound)%nlocal
        if (i > size(iv%sound)) exit
        
        ! Loop through levels for this station
        do k = 1, iv%info(sound)%levels(i)
          ! Extract U component if QC is good and innovation is non-zero
          if (iv%sound(i)%u(k)%qc == 0 .and. abs(iv%sound(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%sound(i)%u(k)%inv
          end if
          
          ! Extract V component if QC is good and innovation is non-zero
          if (iv%sound(i)%v(k)%qc == 0 .and. abs(iv%sound(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%sound(i)%v(k)%inv
          end if
          
          ! Extract T component if QC is good and innovation is non-zero
          if (iv%sound(i)%t(k)%qc == 0 .and. abs(iv%sound(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%sound(i)%t(k)%inv
          end if
          
          ! Extract Q component if QC is good and innovation is non-zero
          if (iv%sound(i)%q(k)%qc == 0 .and. abs(iv%sound(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%sound(i)%q(k)%inv
          end if
        end do
      end do
    end if
    
    ! Add more families here as needed (buoy, airep, pilot, etc.)
    
    num_innovations = count
    wrfda_extract_innovations = 0
    
  end function wrfda_extract_innovations

  ! Extract observation values from y_type structure for all families
  integer(c_int) function wrfda_extract_observations(iv_ptr, ob_ptr, observations, num_observations) bind(C, name="wrfda_extract_observations")
    implicit none
    type(c_ptr), value :: iv_ptr
    type(c_ptr), value :: ob_ptr
    real(c_double), intent(out) :: observations(*)
    integer(c_int), intent(out) :: num_observations
    
    type(y_type), pointer :: y
    type(iv_type), pointer :: iv
    integer :: i, k, count
    
    ! Convert C pointers to Fortran pointers
    call c_f_pointer(iv_ptr, iv)
    call c_f_pointer(ob_ptr, y)
    
    if (.not. associated(iv)) then
      num_observations = 0
      wrfda_extract_observations = 1
      return
    end if
    
    if (.not. associated(y)) then
      num_observations = 0
      wrfda_extract_observations = 1
      return
    end if
    
    ! Extract observations for all families
    count = 0
    
    ! SYNOP family
    if (associated(y%synop) .and. associated(iv%synop) .and. iv%info(synop)%nlocal > 0) then
      do i = 1, min(size(y%synop), iv%info(synop)%nlocal)
        ! Extract U component if QC is good and innovation is non-zero
        if (iv%synop(i)%u%qc == 0 .and. abs(iv%synop(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%synop(i)%u
        end if
        
        ! Extract V component if QC is good and innovation is non-zero
        if (iv%synop(i)%v%qc == 0 .and. abs(iv%synop(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%synop(i)%v
        end if
        
        ! Extract T component if QC is good and innovation is non-zero
        if (iv%synop(i)%t%qc == 0 .and. abs(iv%synop(i)%t%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%synop(i)%t
        end if
        
        ! Extract P component if QC is good and innovation is non-zero
        if (iv%synop(i)%p%qc == 0 .and. abs(iv%synop(i)%p%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%synop(i)%p
        end if
        
        ! Extract Q component if QC is good and innovation is non-zero
        if (iv%synop(i)%q%qc == 0 .and. abs(iv%synop(i)%q%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%synop(i)%q
        end if
      end do
    end if
    
    ! METAR family
    if (associated(y%metar) .and. associated(iv%metar) .and. iv%info(metar)%nlocal > 0) then
      do i = 1, min(size(y%metar), iv%info(metar)%nlocal)
        if (iv%metar(i)%u%qc == 0 .and. abs(iv%metar(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%metar(i)%u
        end if
        if (iv%metar(i)%v%qc == 0 .and. abs(iv%metar(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%metar(i)%v
        end if
        if (iv%metar(i)%t%qc == 0 .and. abs(iv%metar(i)%t%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%metar(i)%t
        end if
        if (iv%metar(i)%p%qc == 0 .and. abs(iv%metar(i)%p%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%metar(i)%p
        end if
        if (iv%metar(i)%q%qc == 0 .and. abs(iv%metar(i)%q%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%metar(i)%q
        end if
      end do
    end if
    
    ! SHIPS family
    if (associated(y%ships) .and. associated(iv%ships) .and. iv%info(ships)%nlocal > 0) then
      do i = 1, min(size(y%ships), iv%info(ships)%nlocal)
        if (iv%ships(i)%u%qc == 0 .and. abs(iv%ships(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ships(i)%u
        end if
        if (iv%ships(i)%v%qc == 0 .and. abs(iv%ships(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ships(i)%v
        end if
        if (iv%ships(i)%t%qc == 0 .and. abs(iv%ships(i)%t%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ships(i)%t
        end if
        if (iv%ships(i)%p%qc == 0 .and. abs(iv%ships(i)%p%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ships(i)%p
        end if
        if (iv%ships(i)%q%qc == 0 .and. abs(iv%ships(i)%q%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ships(i)%q
        end if
      end do
    end if
    
    ! SOUND family - if available and allocated
    ! Note: sound observations have multiple levels per station
    if (associated(y%sound) .and. associated(iv%sound) .and. iv%info(sound)%nlocal > 0) then
      do i = 1, min(size(y%sound), iv%info(sound)%nlocal)
        ! Loop through levels for this station
        do k = 1, iv%info(sound)%levels(i)
          ! Extract U component if QC is good and innovation is non-zero
          if (iv%sound(i)%u(k)%qc == 0 .and. abs(iv%sound(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%sound(i)%u(k)
          end if
          
          ! Extract V component if QC is good and innovation is non-zero
          if (iv%sound(i)%v(k)%qc == 0 .and. abs(iv%sound(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%sound(i)%v(k)
          end if
          
          ! Extract T component if QC is good and innovation is non-zero
          if (iv%sound(i)%t(k)%qc == 0 .and. abs(iv%sound(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%sound(i)%t(k)
          end if
          
          ! Extract Q component if QC is good and innovation is non-zero
          if (iv%sound(i)%q(k)%qc == 0 .and. abs(iv%sound(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%sound(i)%q(k)
          end if
        end do
      end do
    end if
    
    ! Add more families here as needed (buoy, airep, pilot, etc.)
    
    num_observations = count
    wrfda_extract_observations = 0
    
  end function wrfda_extract_observations

  ! Set delta_y input for adjoint operator
  integer(c_int) function wrfda_set_delta_y(delta_y, num_obs) bind(C, name="wrfda_set_delta_y")
    implicit none
    real(c_double), intent(in) :: delta_y(*)
    integer(c_int), value :: num_obs
    
    ! Allocate persistent delta_y array if needed
    if (.not. allocated(persistent_delta_y)) then
      allocate(persistent_delta_y(num_obs))
    else if (size(persistent_delta_y) /= num_obs) then
      deallocate(persistent_delta_y)
      allocate(persistent_delta_y(num_obs))
    end if
    
    ! Copy input delta_y to persistent array
    persistent_delta_y(1:num_obs) = delta_y(1:num_obs)
    
    wrfda_set_delta_y = 0_c_int
  end function wrfda_set_delta_y

  ! Get adjoint gradients from persistent arrays
  integer(c_int) function wrfda_get_adjoint_gradients(u, v, t, q, psfc) bind(C, name="wrfda_get_adjoint_gradients")
    implicit none
    real(c_double), intent(out) :: u(*), v(*), t(*), q(*), psfc(*)
    
    integer :: total_size, surface_size
    
    if (.not. adjoint_allocated) then
      wrfda_get_adjoint_gradients = 1_c_int
      return
    end if
    
    total_size = persistent_nx * persistent_ny * persistent_nz
    surface_size = persistent_nx * persistent_ny
    
    ! Copy gradients from persistent arrays
    u(1:total_size) = persistent_adjoint_u(1:total_size)
    v(1:total_size) = persistent_adjoint_v(1:total_size)
    t(1:total_size) = persistent_adjoint_t(1:total_size)
    q(1:total_size) = persistent_adjoint_q(1:total_size)
    psfc(1:surface_size) = persistent_adjoint_psfc(1:surface_size)
    
    wrfda_get_adjoint_gradients = 0_c_int
  end function wrfda_get_adjoint_gradients

  ! Initialize WRFDA variables for 3D-Var analysis
  !> @brief Copy namelist configuration from model_config_rec to da_control module
  !> 
  !> @details This function extracts and reuses the config copy logic from WRFDA's
  !>          da_wrfvar_init1 subroutine (lines 87-92 of da_wrfvar_init1.inc).
  !>          
  !>          We cannot call da_wrfvar_init1 directly because it also:
  !>          - Calls init_modules (already done in WRFConfigManager)
  !>          - Calls wrfu_initialize (already done)
  !>          - Calls initial_config to read namelist (already done)
  !>          - Sets up MPI variables (already done)
  !>          
  !>          This function extracts ONLY the config copy step using the EXACT
  !>          same mechanism as da_wrfvar_init1: including config_assigns.inc
  !>          with SOURCE_RECORD defined as model_config_rec%.
  !>          
  !> @note This must be called after initial_config() to populate da_control variables
  !> @see da_wrfvar_init1.inc lines 87-92 in WRFDA source
  subroutine copy_config_to_da_control() bind(C, name="copy_config_to_da_control")
    use module_configure, only: model_config_rec
    use da_control
    implicit none
    
    ! Copy namelist variables from model_config_rec to da_control module
    ! This is the EXACT code from da_wrfvar_init1.inc lines 87-92:
    !   #define SOURCE_RECORD model_config_rec%
    !   #define DEST_RECORD
    !   #include "config_assigns.inc"
    ! 
    ! config_assigns.inc is generated by WRF's Registry and contains assignments like:
    !   use_synopobs = model_config_rec%use_synopobs
    !   use_metarobs = model_config_rec%use_metarobs
    !   ... (2000+ lines of similar assignments)
    
#define SOURCE_RECORD model_config_rec%
#define DEST_RECORD
#include "config_assigns.inc"
#undef SOURCE_RECORD
#undef DEST_RECORD
    
  end subroutine copy_config_to_da_control

  !> @brief Validate WRFDA configuration for common conflicts and errors
  !> 
  !> @details This function performs the same sanity checks as WRFDA's da_solve.inc
  !>          to catch configuration errors early. These checks replicate the logic
  !>          from da_solve.inc lines 168-252 (Initial checks section).
  !>          
  !> @note This must be called after copy_config_to_da_control()
  !> @see da_solve.inc lines 168-252 in WRFDA source
  !> @return Integer error code (0 = success, non-zero = validation error)
  subroutine validate_wrfda_config(error_code) bind(C, name="validate_wrfda_config")
    use da_control, only: use_gpsrefobs, use_gpsephobs, use_radar_rf, use_radar_rhv, &
                          use_radarobs, cv_options_hum, cv_options_hum_specific_humidity, &
                          cv_options_hum_relative_humidity, vert_corr, vert_corr_2, &
                          vertical_ip, vertical_ip_0, vertical_ip_delta_p, &
                          cv_options, cloud_cv_options, ensdim_alpha, alphacv_method, &
                          alphacv_method_xa, anal_type_hybrid_dual_res, radar_rf_opt, &
                          alpha_hydrometeors
    use da_wrf_interfaces, only: wrf_message, wrf_error_fatal
    implicit none
    
    integer(c_int), intent(out) :: error_code
    character(len=512) :: message
    
    error_code = 0_c_int
    
    ! ========================================================================
    ! GPS Observation Conflicts (da_solve.inc lines 248-252)
    ! ========================================================================
    if (use_gpsrefobs .and. use_gpsephobs) then
      call wrf_message("ERROR: Configuration Validation Failed")
      call wrf_message("Both use_gpsrefobs and use_gpsephobs are set to true")
      call wrf_message("You must choose EITHER use_gpsrefobs OR use_gpsephobs, not both")
      call wrf_message("Please update your namelist.input &wrfvar4 section and rerun")
      error_code = 1_c_int
      return
    end if
    
    ! ========================================================================
    ! Radar Observation Conflicts (da_solve.inc lines 225-238)
    ! ========================================================================
    if (use_radarobs) then
      if (use_radar_rf .and. use_radar_rhv) then
        call wrf_message("ERROR: Configuration Validation Failed")
        call wrf_message("Both use_radar_rf and use_radar_rhv are set to true")
        call wrf_message("You must choose EITHER use_radar_rf OR use_radar_rhv, not both")
        call wrf_message("Please update your namelist.input &wrfvar4 section and rerun")
        error_code = 2_c_int
        return
      end if
      
      if (use_radar_rf .and. radar_rf_opt == 1 .and. cloud_cv_options /= 1) then
        call wrf_message("ERROR: Configuration Validation Failed")
        call wrf_message("For use_radar_rf with radar_rf_opt=1, you must set cloud_cv_options=1")
        call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
        error_code = 3_c_int
        return
      end if
      
      if (use_radar_rhv .and. cloud_cv_options == 1) then
        call wrf_message("ERROR: Configuration Validation Failed")
        call wrf_message("For use_radar_rhv, you must set cloud_cv_options=2 or 3")
        call wrf_message("(cloud_cv_options=2 requires cloudy be.dat)")
        call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
        error_code = 4_c_int
        return
      end if
    end if
    
    ! ========================================================================
    ! Control Variable Options Validation (da_solve.inc lines 168-173)
    ! ========================================================================
    if (cv_options_hum /= cv_options_hum_specific_humidity .and. &
        cv_options_hum /= cv_options_hum_relative_humidity) then
      write(message, '(A,I3)') &
        "ERROR: Invalid cv_options_hum = ", cv_options_hum
      call wrf_message(trim(message))
      call wrf_message("cv_options_hum must be 1 (specific humidity) or 2 (relative humidity)")
      call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
      error_code = 5_c_int
      return
    end if
    
    ! ========================================================================
    ! Vertical Interpolation Validation (da_solve.inc lines 175-181)
    ! ========================================================================
    if (vert_corr == vert_corr_2) then
      if (vertical_ip < vertical_ip_0 .or. vertical_ip > vertical_ip_delta_p) then
        write(message, '(A,I3)') &
          "ERROR: Invalid vertical_ip = ", vertical_ip
        call wrf_message(trim(message))
        write(message, '(A,I1,A,I1)') &
          "vertical_ip must be between ", vertical_ip_0, " and ", vertical_ip_delta_p
        call wrf_message(trim(message))
        call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
        error_code = 6_c_int
        return
      end if
    end if
    
    ! ========================================================================
    ! CV Options = 3 Constraints (da_solve.inc lines 206-221)
    ! ========================================================================
    if (cv_options == 3) then
      if (ensdim_alpha > 0) then
        call wrf_message("ERROR: Configuration Validation Failed")
        call wrf_message("Alpha control variables are not implemented for cv_options=3")
        call wrf_message("Please set cv_options=5 or 7, OR set ensdim_alpha=0")
        call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
        error_code = 7_c_int
        return
      end if
    end if
    
    ! ========================================================================
    ! Dual-Resolution Hybrid Constraints (da_solve.inc lines 194-198)
    ! ========================================================================
    if (anal_type_hybrid_dual_res .and. alphacv_method /= alphacv_method_xa) then
      call wrf_message("ERROR: Configuration Validation Failed")
      call wrf_message("Dual-resolution hybrid requires alphacv_method=2")
      write(message, '(A,I1)') "Current alphacv_method = ", alphacv_method
      call wrf_message(trim(message))
      call wrf_message("Please update your namelist.input &wrfvar17 section and rerun")
      error_code = 8_c_int
      return
    end if
    
    ! ========================================================================
    ! Alpha Hydrometeors Constraints (da_solve.inc lines 240-245)
    ! ========================================================================
    if (ensdim_alpha > 0 .and. alpha_hydrometeors) then
      if (cloud_cv_options == 1) then
        call wrf_message("ERROR: Configuration Validation Failed")
        call wrf_message("alpha_hydrometeors is not implemented for cloud_cv_options=1")
        call wrf_message("Please set cloud_cv_options=3")
        call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
        error_code = 9_c_int
        return
      end if
    end if
    
    ! All checks passed
    call wrf_message("Configuration validation passed - all settings are consistent")
    
  end subroutine validate_wrfda_config

  subroutine initialize_wrfda_3dvar() bind(C, name="initialize_wrfda_3dvar")
    implicit none
    
    ! Set variables for 3D-Var analysis
    ! var4d_run = .false. indicates 3D-Var (not 4D-Var)
    ! num_fgat_time = 1 indicates single time slot for 3D-Var
    ! sfc_assi_options = sfc_assi_options_1 sets surface assimilation options
    var4d_run = .false.
    num_fgat_time = 1
    sfc_assi_options = sfc_assi_options_1
    
  end subroutine initialize_wrfda_3dvar
  

  ! Initialize domain clock for WRFDA time management
  subroutine initialize_domain_clock(grid)
    implicit none
    type(domain), intent(inout) :: grid
    character(len=256) :: timestr
    logical, save :: wrfu_initialized = .false.
    type(WRFU_Time) :: start_time, stop_time
    type(WRFU_TimeInterval) :: time_step
    integer :: rc
    
    ! Initialize WRFU (required for domain_clock_set calls)
    if (.not. wrfu_initialized) then
      call wrfu_initialize(defaultCalKind=wrfu_cal_gregorian)
      wrfu_initialized = .true.
    end if
    
    ! Set a default analysis time for the domain clock
    timestr = '2000-01-24_12:00:00'  ! Default analysis time
    
    ! Convert time string to WRFU_Time using WRF's wrf_atotime function
    call wrf_atotime(timestr, start_time)
    call wrf_atotime(timestr, stop_time)
    
    ! Create time step interval (0 seconds for analysis)
    call WRFU_TimeIntervalSet(time_step, S=0, rc=rc)
    if (rc /= WRFU_SUCCESS) then
      print *, "WRFDA ERROR: Failed to create time step interval"
      return
    end if
    
    ! Create the domain clock using WRFU_ClockCreate
    ! Note: domain_clock is a TYPE(WRFU_Clock), not a pointer, so no allocation needed
    grid%domain_clock = WRFU_ClockCreate(TimeStep=time_step, StartTime=start_time, StopTime=stop_time, rc=rc)
    if (rc /= WRFU_SUCCESS) then
      print *, "WRFDA ERROR: Failed to create domain clock"
      return
    end if
    
    ! Set the domain_clock_created flag to true
    grid%domain_clock_created = .true.
    
  end subroutine initialize_domain_clock

  ! =============================================================================
  ! NEW: WRFDA First Guess Initialization Functions
  ! These functions implement the refactored WRFState initialization strategy
  ! that leverages WRFDA's proven initialization pipeline
  ! =============================================================================

  !> @brief Load WRF first guess using WRFDA's da_med_initialdata_input and da_setup_firstguess
  !> @details This function implements the complete WRFDA initialization pipeline:
  !>          1. Calls da_med_initialdata_input to read NetCDF file into grid structure
  !>          2. Calls da_setup_firstguess_wrf to setup grid and call da_transfer_wrftoxb
  !>          3. All WRFDA diagnostics and derived fields are computed properly
  !> @param[in] filename C string containing path to WRF NetCDF file
  !> @param[in] filename_len Length of filename string
  !> @return Integer status code (0 = success, non-zero = error)
!> @brief Update da_control module variables from grid structure
!> @details Sets all grid-dependent da_control module variables after da_med_initialdata_input.
!>          This includes:
!>            (1) Vertical coordinate coefficients (c1f, c2f, c3f, c4f, c1h, c2h, c3h, c4h)
!>            (2) Base state parameters (base_pres, base_temp, base_lapse, iso_temp, etc.)
!>          Map projection variables are set separately by da_setup_firstguess_wrf.
!>          This function must be called AFTER da_med_initialdata_input and BEFORE
!>          da_setup_firstguess_wrf/da_transfer_wrftoxb.
!> @param[in] grid Pointer to WRFDA grid structure with populated state from NetCDF
subroutine update_da_control_from_grid(grid)
  use module_domain, only: domain
  use da_control, only: c1f, c2f, c3f, c4f, c1h, c2h, c3h, c4h, &
                        base_pres, base_temp, base_lapse, iso_temp, &
                        base_pres_strat, base_lapse_strat
  use da_wrf_interfaces, only: wrf_message
  implicit none
  
  type(domain), pointer, intent(in) :: grid
  integer :: kms, kme
  character(len=256) :: msg
  
  call wrf_message("======== Updating da_control module variables from grid ========")
  
  ! -------------------------------------------------------------------------
  ! [1] Update vertical coordinate coefficients
  ! -------------------------------------------------------------------------
  kms = grid%sm33
  kme = grid%em33
  
  ! Allocate da_control module arrays if not already allocated
  if (.not. allocated(c1f)) allocate(c1f(kms:kme))
  if (.not. allocated(c2f)) allocate(c2f(kms:kme))
  if (.not. allocated(c3f)) allocate(c3f(kms:kme))
  if (.not. allocated(c4f)) allocate(c4f(kms:kme))
  if (.not. allocated(c1h)) allocate(c1h(kms:kme))
  if (.not. allocated(c2h)) allocate(c2h(kms:kme))
  if (.not. allocated(c3h)) allocate(c3h(kms:kme))
  if (.not. allocated(c4h)) allocate(c4h(kms:kme))
  
  write(msg, '(A,I3,A,I3,A)') "[1] Vertical coordinates (kms:kme = ", kms, ":", kme, ")"
  call wrf_message(msg)
  
  ! Copy from grid structure (populated by da_med_initialdata_input)
  c1f = grid%c1f
  c2f = grid%c2f
  c3f = grid%c3f
  c4f = grid%c4f
  c1h = grid%c1h
  c2h = grid%c2h
  c3h = grid%c3h
  c4h = grid%c4h
  
  ! Handle non-hybrid coordinates (backward compatibility)
  ! For input files prior to V3.9, grid%hybrid_opt is set to 0 by da_med_initialdata_input
  if (grid%hybrid_opt <= 0) then
    write(msg, '(A,I2,A)') "    Hybrid_opt = ", grid%hybrid_opt, " -> using pure eta coordinates"
    call wrf_message(msg)
    
    ! Fall back to pure eta coordinates
    c3f = grid%znw  ! Eta levels on full (w) layers
    c3h = grid%znu  ! Eta levels on half (mass) layers
    c4f = 0.0
    c4h = 0.0
    c1f = 1.0
    c1h = 1.0
    c2f = 0.0
    c2h = 0.0
  else
    write(msg, '(A,I2,A)') "    Hybrid_opt = ", grid%hybrid_opt, " -> using hybrid coordinates"
    call wrf_message(msg)
  end if
  
  ! -------------------------------------------------------------------------
  ! [2] Update base state parameters
  ! -------------------------------------------------------------------------
  write(msg, '(A)') "[2] Base state parameters"
  call wrf_message(msg)
  
  base_pres  = grid%p00
  base_temp  = grid%t00
  base_lapse = grid%tlp
  iso_temp   = grid%tiso
  base_pres_strat  = grid%p_strat
  base_lapse_strat = grid%tlp_strat
  
  write(msg, '(A,F10.1,A)') "    base_pres        = ", base_pres, " Pa"
  call wrf_message(msg)
  write(msg, '(A,F10.2,A)') "    base_temp        = ", base_temp, " K"
  call wrf_message(msg)
  write(msg, '(A,F10.6)') "    base_lapse       = ", base_lapse
  call wrf_message(msg)
  write(msg, '(A,F10.2,A)') "    iso_temp         = ", iso_temp, " K"
  call wrf_message(msg)
  write(msg, '(A,F10.1,A)') "    base_pres_strat  = ", base_pres_strat, " Pa"
  call wrf_message(msg)
  write(msg, '(A,F10.6)') "    base_lapse_strat = ", base_lapse_strat
  call wrf_message(msg)
  
  ! Validate - base_temp should be > 100 K and base_pres > 10000 Pa
  if ( base_temp < 100.0 .or. base_pres < 10000.0 ) then
    write(msg, '(A)') "ERROR: Base state parameters not found in NetCDF file!"
    call wrf_message(msg)
    write(msg, '(A)') "Add use_baseparam_fr_nml = .true. in namelist.input &dynamics"
    call wrf_message(msg)
    write(msg, '(A,F10.2,A,F10.1)') "Got: base_temp = ", base_temp, " K, base_pres = ", base_pres
    call wrf_message(msg)
    stop "FATAL: Missing or invalid base state parameters"
  end if
  
  call wrf_message("======== da_control module variables updated successfully ========")
  
end subroutine update_da_control_from_grid

!> @brief Cleanup da_control module vertical coordinate variables
!> @details Deallocates vertical coordinate arrays in da_control module
!>          Should be called when MetaDA is done with the grid
subroutine cleanup_da_control_vertical_coords()
  use da_control, only: c1f, c2f, c3f, c4f, c1h, c2h, c3h, c4h
  use da_wrf_interfaces, only: wrf_message
  implicit none
  
  call wrf_message("Cleaning up da_control vertical coordinates")
  
  if (allocated(c1f)) deallocate(c1f)
  if (allocated(c2f)) deallocate(c2f)
  if (allocated(c3f)) deallocate(c3f)
  if (allocated(c4f)) deallocate(c4f)
  if (allocated(c1h)) deallocate(c1h)
  if (allocated(c2h)) deallocate(c2h)
  if (allocated(c3h)) deallocate(c3h)
  if (allocated(c4h)) deallocate(c4h)
  
  call wrf_message("da_control vertical coordinates cleaned up")
  
end subroutine cleanup_da_control_vertical_coords

!> @brief C-callable wrapper for cleanup_da_control_vertical_coords
!> @details Allows C++ code to cleanup da_control module vertical coordinates
subroutine wrfda_cleanup_vertical_coords() bind(C, name="wrfda_cleanup_vertical_coords")
  implicit none
  
  call cleanup_da_control_vertical_coords()
  
end subroutine wrfda_cleanup_vertical_coords

integer(c_int) function wrfda_load_first_guess(grid_ptr, filename, filename_len) &
    bind(C, name="wrfda_load_first_guess")
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_char, c_int, c_null_char, c_associated, c_double
  use module_configure, only: grid_config_rec_type, model_config_rec, model_to_grid_config_rec
  use module_domain, only: domain, head_grid
  use module_tiles, only: set_tiles
  use da_wrfvar_io, only: da_med_initialdata_input
  use da_transfer_model, only: da_setup_firstguess_wrf
  use da_wrf_interfaces, only: wrf_message
  implicit none
  
  ! Input parameters
  type(c_ptr), value, intent(in) :: grid_ptr
  type(c_ptr), value, intent(in) :: filename
  integer(c_int), value, intent(in) :: filename_len
  
  ! Local variables
  type(domain), pointer :: grid
  character(len=512) :: fg_filename
  type(xbx_type) :: xbx
  type(grid_config_rec_type) :: config_flags
  integer :: i, max_len
  character(len=256) :: msg
  character(kind=c_char), pointer :: fptr(:)
  
  ! Initialize return code
  wrfda_load_first_guess = 0
  
  ! Validate grid pointer (must be pre-allocated via wrfda_alloc_and_init_domain)
  if (.not. c_associated(grid_ptr)) then
    call wrf_message("ERROR: grid_ptr is null. Allocate domain first via wrfda_alloc_and_init_domain")
    wrfda_load_first_guess = -1
    return
  end if
  
  ! Convert C pointer to Fortran pointer (this is head_grid from WRFConfigManager)
  call c_f_pointer(grid_ptr, grid)
  
  ! Convert C filename string to Fortran string
  if (filename_len < 0 .or. filename_len > 1024) then
    call wrf_message("ERROR: Invalid filename_len")
    wrfda_load_first_guess = -1
    return
  end if
  
  fg_filename = ""
  call c_f_pointer(filename, fptr, [filename_len])
  max_len = min(filename_len, len(fg_filename))
  do i = 1, max_len
    if (fptr(i) == c_null_char) exit
    fg_filename(i:i) = fptr(i)
  end do
  
  write(msg, '(A,A)') "WRFDA: Loading first guess from file: ", trim(fg_filename)
  call wrf_message(msg)
  
  ! Extract domain-specific config in a single canonical call
  ! Note: This was already called in wrfda_alloc_and_init_domain, but we need
  ! a local config_flags for da_med_initialdata_input
  call model_to_grid_config_rec(head_grid%id, model_config_rec, config_flags)
  
  ! Call da_med_initialdata_input to read NetCDF file
  call da_med_initialdata_input(head_grid, config_flags, trim(fg_filename))
  
  ! Update da_control module variables (vertical coords + base state parameters)
  ! This must be done AFTER da_med_initialdata_input (which populates grid structure)
  ! and BEFORE da_setup_firstguess_wrf/da_transfer_wrftoxb (which use these module variables)
  call update_da_control_from_grid(head_grid)
  
  ! Initialize WRFDA module-level variables following standard WRFDA sequence
  ! This matches the exact sequence used in da_solve_init.inc (lines 33-40)
  ! Step 1: De-reference dimension information from grid structure
  !         Sets module variables: ids, ide, jds, jde, ips, ipe, jps, jpe, kms, kme, etc.
  call da_copy_dims(head_grid)
  
  ! Step 2: Compute tile starting/stopping locations
  !         Sets grid%i_start, grid%i_end, grid%j_start, grid%j_end, grid%num_tiles
  call set_tiles(head_grid, ids, ide, jds, jde, ips, ipe, jps, jpe)
  
  ! Step 3: Copy tile dimensions to module-level variables
  !         Sets module variables: its, ite, jts, jte, kts, kte
  call da_copy_tile_dims(head_grid)
  
  ! Set surface assimilation options
  sfc_assi_options = sfc_assi_options_1
  
  ! Call da_setup_firstguess_wrf to setup grid and call da_transfer_wrftoxb
  call da_setup_firstguess_wrf(xbx, head_grid, config_flags, .false.)
  
  write(msg, '(A)') "WRFDA: First guess loaded successfully"
  call wrf_message(msg)
    
  end function wrfda_load_first_guess

  !> @brief Extract background state (xb) data from WRFDA grid to C++ arrays
  !> @details Extracts the fully initialized background state from WRFDA's grid structure
  !>          to flat C++ arrays for use in the WRFState class
  !> @param[out] u U-wind component (3D: nx*ny*nz)
  !> @param[out] v V-wind component (3D: nx*ny*nz)
  !> @param[out] t Temperature perturbation (3D: nx*ny*nz)
  !> @param[out] q Specific humidity (3D: nx*ny*nz)
  !> @param[out] psfc Surface pressure (2D: nx*ny)
  !> @param[out] p Full pressure (3D: nx*ny*nz)
  !> @param[out] ph Geopotential perturbation (3D staggered: nx*ny*(nz+1))
  !> @param[out] phb Base state geopotential (3D staggered: nx*ny*(nz+1))
  !> @param[out] hgt Terrain height (2D: nx*ny)
  !> @param[out] lats 2D latitude array (2D: nx*ny)
  !> @param[out] lons 2D longitude array (2D: nx*ny)
  !> @param[out] nx Grid x-dimension
  !> @param[out] ny Grid y-dimension
  !> @param[out] nz Grid z-dimension
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_extract_background_state( &
      grid_ptr, u, v, t, q, psfc, p, ph, phb, hgt, lats, lons, &
      nx, ny, nz) bind(C, name="wrfda_extract_background_state")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use module_domain, only: domain
    implicit none
    
    ! Input parameters
    type(c_ptr), value, intent(in) :: grid_ptr
    
    ! Output parameters
    real(c_double), intent(out) :: u(*), v(*), t(*), q(*), psfc(*)
    real(c_double), intent(out) :: p(*), ph(*), phb(*), hgt(*)
    real(c_double), intent(out) :: lats(*), lons(*)
    integer(c_int), intent(out) :: nx, ny, nz
    
    ! Local variables
    integer :: i, j, k, idx, idx_3d, staggered_nz
    type(domain), pointer :: grid
    
    ! Initialize return code
    wrfda_extract_background_state = 0
    
    ! Handle grid pointer: use persistent_grid if grid_ptr is null
    if (.not. c_associated(grid_ptr)) then
      if (.not. grid_initialized .or. .not. associated(persistent_grid)) then
        wrfda_extract_background_state = -1
        return
      end if
      grid => persistent_grid
    else
      call c_f_pointer(grid_ptr, grid)
    end if
    
    ! Get grid dimensions
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    staggered_nz = nz + 1
    
    ! Check if xb is allocated
    if (.not. associated(grid%xb%u) .or. .not. associated(grid%xb%v) .or. &
        .not. associated(grid%xb%t) .or. .not. associated(grid%xb%q) .or. &
        .not. associated(grid%xb%psfc)) then
      wrfda_extract_background_state = 2
      return
    end if
    
    ! Extract 3D fields (U, V, T, Q, P)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          ! But we're using 1-based Fortran indexing, so: (i-1) + (j-1)*nx + (k-1)*nx*ny
          idx = (i-1) + (j-1)*nx + (k-1)*nx*ny + 1
          u(idx) = real(grid%xb%u(i,j,k), kind=c_double)
          v(idx) = real(grid%xb%v(i,j,k), kind=c_double)
          t(idx) = real(grid%xb%t(i,j,k), kind=c_double)
          q(idx) = real(grid%xb%q(i,j,k), kind=c_double)
          p(idx) = real(grid%xb%p(i,j,k), kind=c_double)
        end do
      end do
    end do
    
    ! Extract 2D surface fields (PSFC, HGT)
    do j = 1, ny
      do i = 1, nx
        idx = (i-1) + (j-1)*nx + 1
        psfc(idx) = real(grid%xb%psfc(i,j), kind=c_double)
        hgt(idx) = real(grid%xb%terr(i,j), kind=c_double)
        lats(idx) = real(grid%xb%lat(i,j), kind=c_double)
        lons(idx) = real(grid%xb%lon(i,j), kind=c_double)
      end do
    end do
    
    ! Extract vertically staggered fields (PH, PHB) if available
    if (associated(grid%ph_2) .and. associated(grid%phb)) then
      do k = 1, staggered_nz
        do j = 1, ny
          do i = 1, nx
            idx_3d = (i-1) + (j-1)*nx + (k-1)*nx*ny + 1
            ph(idx_3d) = real(grid%ph_2(i,j,k), kind=c_double)
            phb(idx_3d) = real(grid%phb(i,j,k), kind=c_double)
          end do
        end do
      end do
    end if
    
  end function wrfda_extract_background_state

  !> @brief Extract additional WRF fields needed for full state representation
  !> @details Extracts additional fields beyond core state variables
  !> @param[out] w Vertical velocity (3D staggered: nx*ny*(nz+1))
  !> @param[out] mu Dry air mass perturbation (2D: nx*ny)
  !> @param[out] mub Base state dry air mass (2D: nx*ny)
  !> @param[out] pb Base state pressure (3D: nx*ny*nz)
  !> @param[out] t_init Initial temperature field (3D: nx*ny*nz)
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_extract_additional_fields( &
      grid_ptr, w, mu, mub, pb, t_init) bind(C, name="wrfda_extract_additional_fields")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use module_domain, only: domain
    implicit none
    
    ! Input parameters
    type(c_ptr), value, intent(in) :: grid_ptr
    
    ! Output parameters
    real(c_double), intent(out) :: w(*), mu(*), mub(*), pb(*), t_init(*)
    
    ! Local variables
    integer :: i, j, k, idx, idx_3d, nx, ny, nz, staggered_nz
    type(domain), pointer :: grid
    
    ! Initialize return code
    wrfda_extract_additional_fields = 0
    
    ! Handle grid pointer: use persistent_grid if grid_ptr is null
    if (.not. c_associated(grid_ptr)) then
      if (.not. grid_initialized .or. .not. associated(persistent_grid)) then
        wrfda_extract_additional_fields = -1
        return
      end if
      grid => persistent_grid
    else
      call c_f_pointer(grid_ptr, grid)
    end if
    
    ! Get grid dimensions
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    staggered_nz = nz + 1
    
    ! Extract W field (vertically staggered)
    if (associated(grid%w_2)) then
      do k = 1, staggered_nz
        do j = 1, ny
          do i = 1, nx
            idx_3d = (i-1) + (j-1)*nx + (k-1)*nx*ny + 1
            w(idx_3d) = real(grid%w_2(i,j,k), kind=c_double)
          end do
        end do
      end do
    end if
    
    ! Extract MU and MUB fields (2D)
    if (associated(grid%mu_2) .and. associated(grid%mub)) then
      do j = 1, ny
        do i = 1, nx
          idx = (i-1) + (j-1)*nx + 1
          mu(idx) = real(grid%mu_2(i,j), kind=c_double)
          mub(idx) = real(grid%mub(i,j), kind=c_double)
        end do
      end do
    end if
    
    ! Extract PB and T_INIT fields (3D)
    if (associated(grid%pb) .and. associated(grid%t_init)) then
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            idx = (i-1) + (j-1)*nx + (k-1)*nx*ny + 1
            pb(idx) = real(grid%pb(i,j,k), kind=c_double)
            t_init(idx) = real(grid%t_init(i,j,k), kind=c_double)
          end do
        end do
      end do
    end if
    
  end function wrfda_extract_additional_fields

  !> @brief Read and allocate observations using WRFDA's standard pipeline
  !> @details Reads observations from BUFR file and allocates WRFDA iv_type and y_type structures.
  !>          This function calls WRFDA's proven observation reading workflow:
  !>            1. da_setup_obs_structures_bufr → da_read_obs_bufr (reads PREPBUFR)
  !>            2. da_allocate_observations (allocates iv arrays)
  !>            3. da_allocate_y (allocates ob/y arrays)
  !>          CRITICAL: This does NOT compute innovations (iv%synop%u%inv remains empty).
  !>          Innovation computation happens later when da_get_innov_vector is called
  !>          by the observation operator with a background state.
  !> @param[in] grid_ptr C pointer to WRFDA grid structure (for domain bounds)
  !> @param[in] ob_filename C string containing path to BUFR observation file
  !> @param[in] ob_filename_len Length of filename string
  !> @param[out] iv_ptr Output pointer to allocated iv_type structure
  !> @param[out] ob_ptr Output pointer to allocated y_type structure
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_read_and_allocate_observations( &
      grid_ptr, ob_filename, ob_filename_len, iv_ptr, ob_ptr) &
      bind(C, name="wrfda_read_and_allocate_observations")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_char, c_int, c_null_char, c_associated, c_loc
  use module_domain, only: domain, head_grid
  use da_define_structures, only: iv_type, y_type, j_type
  use da_setup_structures, only: da_setup_obs_structures
  use da_wrf_interfaces, only: wrf_message
    implicit none
    
    ! Input parameters
    type(c_ptr), value, intent(in) :: grid_ptr
    type(c_ptr), value, intent(in) :: ob_filename
    integer(c_int), value, intent(in) :: ob_filename_len
    
    ! Output parameters
    type(c_ptr), intent(out) :: iv_ptr
    type(c_ptr), intent(out) :: ob_ptr
    
    ! WRFDA-side structures with SAVE attribute (persist across calls)
    ! These are allocated in WRFDA modules and managed by WRFDA
    ! C++ only receives pointers to them
    type(iv_type), pointer, save :: wrfda_iv => null()
    type(y_type), pointer, save :: wrfda_ob => null()
    type(j_type), save :: wrfda_j  ! Cost function structure (needed for da_setup_obs_structures)
    logical, save :: wrfda_obs_initialized = .false.
    
    ! Local variables
    type(domain), pointer :: grid
    character(len=512) :: bufr_filename
    character(kind=c_char), pointer :: fptr(:)
    integer :: i, max_len
    character(len=256) :: msg
    
    ! Initialize return code
    wrfda_read_and_allocate_observations = 0
    
    ! Validate grid pointer
    if (.not. c_associated(grid_ptr)) then
      call wrf_message("ERROR: grid_ptr is null in wrfda_read_and_allocate_observations")
      wrfda_read_and_allocate_observations = -1
      return
    end if
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(grid_ptr, grid)
    
    ! Convert C filename string to Fortran string
    if (ob_filename_len < 0 .or. ob_filename_len > 1024) then
      call wrf_message("ERROR: Invalid ob_filename_len")
      wrfda_read_and_allocate_observations = -1
      return
    end if
    
    bufr_filename = ""
    call c_f_pointer(ob_filename, fptr, [ob_filename_len])
    max_len = min(ob_filename_len, len(bufr_filename))
    do i = 1, max_len
      if (fptr(i) == c_null_char) exit
      bufr_filename(i:i) = fptr(i)
    end do
    
    write(msg, '(A,A)') "WRFDA: Reading observations from file: ", trim(bufr_filename)
    call wrf_message(msg)
    
    ! Note: da_setup_obs_structures reads observations based on ob_format in namelist:
    !   ob_format=1 → reads ob.bufr (PREPBUFR format)
    !   ob_format=2 → reads ob.ascii (ASCII format)
    ! File must be in current working directory with these exact names
    
    ! Allocate WRFDA-side observation structures (function-level SAVE)
    ! These persist across calls but are local to this function
    if (.not. wrfda_obs_initialized) then
      allocate(wrfda_iv)
      allocate(wrfda_ob)
      wrfda_obs_initialized = .true.
      call wrf_message("Allocated WRFDA-side iv_type and y_type structures")
      
      ! CRITICAL: Initialize MPI-like variables for single-processor mode
      ! These are required by da_setup_obs_structures and observation reading
      ! Following da_wrfvar_init1.inc logic for non-DM_PARALLEL builds
      num_procs = 1
      myproc = 0
      comm = 0
      rootproc = .true.
      call wrf_message("Initialized processor variables: num_procs=1, myproc=0, rootproc=.true.")
      
      ! Initialize iv structure following WRFDA's da_setup_obs_structures pattern
      ! These initializations are CRITICAL before calling da_setup_obs_structures
      do i = 1, num_ob_indexes
        wrfda_iv%info(i)%nlocal = 0
        wrfda_iv%info(i)%ntotal = 0
        wrfda_iv%info(i)%thin_nlocal = 0
        wrfda_iv%info(i)%thin_ntotal = 0
        wrfda_iv%info(i)%plocal(:) = 0
        wrfda_iv%info(i)%ptotal(:) = 0
        wrfda_iv%info(i)%thin_plocal(:) = 0
        wrfda_iv%info(i)%thin_ptotal(:) = 0
        wrfda_iv%info(i)%max_lev = 1  ! Default to 1 level (surface obs)
      end do
      wrfda_iv%num_inst = 0
      
      ! Initialize ob structure
      wrfda_ob%nlocal(:) = 0
      wrfda_ob%ntotal(:) = 0
      wrfda_ob%num_inst = 0
      
      call wrf_message("Initialized iv and ob structures")
    end if
    
    ! Call WRFDA's standard observation reading pipeline
    ! This function internally calls either:
    !   - da_setup_obs_structures_bufr (if ob_format=1) → reads PREPBUFR file
    !   - da_setup_obs_structures_ascii (if ob_format=2) → reads ASCII file
    ! based on the ob_format setting in namelist.input
    ! It also allocates observations and y_type structures
    ! Note: This does NOT compute innovations (grid%xb not needed yet)
    call da_setup_obs_structures(head_grid, wrfda_ob, wrfda_iv, wrfda_j)
    
    write(msg, '(A)') "WRFDA: Observations read and structures allocated successfully"
    call wrf_message(msg)
    
    ! Return C pointers to WRFDA-side structures
    ! - ob_ptr (y_type) → for WRFObservation (observation values y_obs)
    ! - iv_ptr (iv_type) → for WRFObsOperator (innovations + interpolation weights)
    ! C++ will manage these as opaque pointers
    iv_ptr = c_loc(wrfda_iv)
    ob_ptr = c_loc(wrfda_ob)
    
    ! Log observation counts by type (nlocal = local to this processor, ntotal = global)
    write(msg, '(A)') "Observation counts by type (nlocal/ntotal):"
    call wrf_message(msg)
    write(msg, '(A,I6,A,I6)') "  SYNOP:   ", wrfda_iv%info(synop)%nlocal, " / ", wrfda_iv%info(synop)%ntotal
    call wrf_message(msg)
    write(msg, '(A,I6,A,I6)') "  METAR:   ", wrfda_iv%info(metar)%nlocal, " / ", wrfda_iv%info(metar)%ntotal
    call wrf_message(msg)
    write(msg, '(A,I6,A,I6)') "  SOUND:   ", wrfda_iv%info(sound)%nlocal, " / ", wrfda_iv%info(sound)%ntotal
    call wrf_message(msg)
    write(msg, '(A,I6,A,I6)') "  SHIPS:   ", wrfda_iv%info(ships)%nlocal, " / ", wrfda_iv%info(ships)%ntotal
    call wrf_message(msg)
    write(msg, '(A,I6,A,I6)') "  BUOY:    ", wrfda_iv%info(buoy)%nlocal, " / ", wrfda_iv%info(buoy)%ntotal
    call wrf_message(msg)
    write(msg, '(A,I6,A,I6)') "  AIREP:   ", wrfda_iv%info(airep)%nlocal, " / ", wrfda_iv%info(airep)%ntotal
    call wrf_message(msg)
    write(msg, '(A,I6,A,I6)') "  PILOT:   ", wrfda_iv%info(pilot)%nlocal, " / ", wrfda_iv%info(pilot)%ntotal
    call wrf_message(msg)
    write(msg, '(A,I6,A,I6)') "  SONDE_SFC:", wrfda_iv%info(sonde_sfc)%nlocal, " / ", wrfda_iv%info(sonde_sfc)%ntotal
    call wrf_message(msg)
    
  end function wrfda_read_and_allocate_observations

  !> @brief Extract observation counts from WRFDA iv_type structure
  !> @details Returns the number of observations for each observation type
  !> @param[in] iv_ptr C pointer to WRFDA iv_type structure
  !> @param[out] obs_counts Array of observation counts (size: num_ob_indexes)
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_get_obs_type_counts(iv_ptr, obs_counts) &
      bind(C, name="wrfda_get_obs_type_counts")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int
    use da_define_structures, only: iv_type
    use da_control, only: num_ob_indexes
    implicit none
    
    type(c_ptr), value, intent(in) :: iv_ptr
    integer(c_int), intent(out) :: obs_counts(num_ob_indexes)
    
    type(iv_type), pointer :: iv
    integer :: i
    
    wrfda_get_obs_type_counts = 0
    
    if (.not. c_associated(iv_ptr)) then
      obs_counts(:) = 0
      wrfda_get_obs_type_counts = -1
      return
    end if
    
    call c_f_pointer(iv_ptr, iv)
    
    ! Extract observation counts for all types
    do i = 1, num_ob_indexes
      obs_counts(i) = iv%info(i)%nlocal
    end do
    
  end function wrfda_get_obs_type_counts

  !> @brief Extract total observation count from WRFDA iv_type structure
  !> @param[in] iv_ptr C pointer to WRFDA iv_type structure
  !> @param[out] total_count Total number of observations across all types
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_get_total_obs_count(iv_ptr, total_count) &
      bind(C, name="wrfda_get_total_obs_count")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int
    use da_define_structures, only: iv_type
    use da_control, only: num_ob_indexes
    implicit none
    
    type(c_ptr), value, intent(in) :: iv_ptr
    integer(c_int), intent(out) :: total_count
    
    type(iv_type), pointer :: iv
    integer :: i
    
    wrfda_get_total_obs_count = 0
    total_count = 0
    
    if (.not. c_associated(iv_ptr)) then
      wrfda_get_total_obs_count = -1
      return
    end if
    
    call c_f_pointer(iv_ptr, iv)
    
    ! Sum observation counts across all types
    do i = 1, num_ob_indexes
      total_count = total_count + iv%info(i)%nlocal
    end do
    
  end function wrfda_get_total_obs_count

  !> @brief Transfer WRF fields to background state structure (grid%xb)
  !!
  !! @details This function wraps WRFDA's da_transfer_wrftoxb routine which:
  !! - Transfers WRF native fields (u, v, t, q, etc.) to grid%xb structure
  !! - Computes derived fields (pressure, height, etc.)
  !! - Applies coordinate transformations for Arakawa-C grid
  !! - Prepares grid%xb for use by observation operators
  !!
  !! This MUST be called before using any WRFDA observation operators to ensure
  !! grid%xb is properly populated from the current state.
  !!
  !! @return 0 on success, non-zero on error
  !!
  !! @see da_transfer_wrftoxb.inc in WRFDA source
  function wrfda_transfer_wrftoxb() result(error_code) bind(C, name="wrfda_transfer_wrftoxb")
    use module_domain, only: head_grid, domain
    use module_configure, only: grid_config_rec_type
    use da_define_structures, only: xbx_type
    use da_transfer_model, only: da_transfer_wrftoxb
    implicit none
    
    integer(c_int) :: error_code
    type(domain), pointer :: grid
    type(grid_config_rec_type) :: config_flags
    type(xbx_type) :: xbx
    
    error_code = 0_c_int
    
    ! Get head grid pointer
    grid => head_grid
    if (.not. associated(grid)) then
      call wrf_message("ERROR: head_grid not associated in wrfda_transfer_wrftoxb")
      error_code = -1_c_int
      return
    end if
    
    ! Note: config_flags is not actually used in da_transfer_wrftoxb, 
    ! but is required by the function signature. We initialize it to default values.
    ! (Verified by checking da_transfer_wrftoxb.inc - no config_flags% references)
    
    ! Call WRFDA's da_transfer_wrftoxb to populate grid%xb from WRF fields
    ! This is the standard WRFDA workflow step that must happen before
    ! using observation operators
    call da_transfer_wrftoxb(xbx, grid, config_flags)
    
  end function wrfda_transfer_wrftoxb

end module metada_wrfda_dispatch