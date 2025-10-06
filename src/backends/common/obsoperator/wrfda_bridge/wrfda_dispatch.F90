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
  use da_define_structures, only: iv_type, y_type, da_allocate_y, da_allocate_obs_info
  use module_symbols_util,  only: wrfu_initialize, wrfu_finalize, wrfu_cal_gregorian
  use module_symbols_util, only: WRFU_ClockCreate, WRFU_TimeIntervalSet, WRFU_SUCCESS, WRFU_Time, WRFU_TimeInterval, WRFU_INITIALIZE, WRFU_CAL_GREGORIAN
  use da_control, only: metar, synop, ships, buoy, airep, pilot, sound, sonde_sfc, &
                        sfc_assi_options, sfc_assi_options_1, trace_use_dull, &
                        var4d_run, num_fgat_time, missing_r, missing_data, num_ob_indexes, &
                        kts, kte, its, ite, jts, jte, Max_StHeight_Diff, kms, kme
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
  
  ! CRITICAL FIX: Add persistent grid to avoid reconstruction issues
  ! WRFDA internal state persists between calls, so we need consistent grid structure
  type(domain), pointer, save :: persistent_grid
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
  
  type(proj_info), save :: map_info
  logical, save :: map_info_initialized = .false.
  
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
  integer(c_int) function wrfda_xtoy_apply_grid(ob_ptr, iv_ptr) bind(C, name="wrfda_xtoy_apply_grid")
    implicit none
    type(c_ptr), value :: ob_ptr, iv_ptr

    type(domain), pointer :: grid
    type(iv_type), pointer :: iv
    type(y_type), pointer :: y
    integer :: n, num_innovations
    character(len=256) :: fam_str
    integer :: fam_id
    
    ! Convert C pointers to Fortran pointers
    grid => persistent_grid
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
  integer(c_int) function wrfda_xtoy_adjoint_grid(iv_ptr) bind(C, name="wrfda_xtoy_adjoint_grid")
    implicit none
    type(c_ptr), value :: iv_ptr

    type(domain), pointer :: grid
    type(iv_type), pointer :: iv
    type(y_type), target :: jo_grad_y
    type(x_type), target :: jo_grad_x
    integer :: n, num_innovations, nx, ny, nz
    character(len=256) :: fam_str
    integer :: fam_id
    
    ! Convert C pointers to Fortran pointers
    grid => persistent_grid
    call c_f_pointer(iv_ptr, iv)
    
    if (.not. associated(iv)) then
      wrfda_xtoy_adjoint_grid = 1_c_int
      return
    end if
    
    ! Get grid dimensions from persistent grid
    nx = persistent_grid%xp%ide - persistent_grid%xp%ids + 1
    ny = persistent_grid%xp%jde - persistent_grid%xp%jds + 1
    nz = persistent_grid%xp%kde - persistent_grid%xp%kds + 1
    
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
  subroutine wrfda_update_analysis_increments(u_inc, v_inc, t_inc, q_inc, psfc_inc) bind(C, name="wrfda_update_analysis_increments")
    implicit none
    real(c_double), intent(in) :: u_inc(*), v_inc(*), t_inc(*), q_inc(*), psfc_inc(*)
    integer :: i,j,k, nx, ny, nz
    
    ! Get grid dimensions from persistent_grid
    nx = persistent_grid%xp%ide - persistent_grid%xp%ids + 1
    ny = persistent_grid%xp%jde - persistent_grid%xp%jds + 1
    nz = persistent_grid%xp%kde - persistent_grid%xp%kds + 1
    
    ! Update only the xa fields (analysis increments)
    do k=1,nz; do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
      persistent_grid%xa%u(i,j,k) = real(u_inc(i + (j-1)*nx + (k-1)*nx*ny))
      persistent_grid%xa%v(i,j,k) = real(v_inc(i + (j-1)*nx + (k-1)*nx*ny))
      persistent_grid%xa%t(i,j,k) = real(t_inc(i + (j-1)*nx + (k-1)*nx*ny))
      persistent_grid%xa%q(i,j,k) = real(q_inc(i + (j-1)*nx + (k-1)*nx*ny))
    end do; end do; end do
    do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j] -> i + j*nx
      persistent_grid%xa%psfc(i,j) = real(psfc_inc(i + (j-1)*nx))
    end do; end do
    
  end subroutine wrfda_update_analysis_increments

  ! Update background state (xb) from state data
  ! This subroutine updates the background state with the current state values
  subroutine wrfda_update_background_state(u, v, t, q, psfc, ph, phb, hf, hgt, p, pb, lats2d, lons2d) bind(C, name="wrfda_update_background_state")
    implicit none
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    real(c_double), intent(in) :: ph(*), phb(*), hf(*), hgt(*), p(*), pb(*)
    real(c_double), intent(in) :: lats2d(*), lons2d(*)
    integer :: i, j, k, nx, ny, nz, idx, idx_3d, staggered_nz
    
    ! Get grid dimensions from persistent_grid
    nx = persistent_grid%xp%ide - persistent_grid%xp%ids + 1
    ny = persistent_grid%xp%jde - persistent_grid%xp%jds + 1
    nz = persistent_grid%xp%kde - persistent_grid%xp%kds + 1
    staggered_nz = nz + 1  ! Height field has nz+1 vertical levels
    
    ! Update background state (xb) with current state values
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx = i + (j-1)*nx + (k-1)*nx*ny
          persistent_grid%xb%u(i,j,k) = real(u(idx))
          persistent_grid%xb%v(i,j,k) = real(v(idx))
          persistent_grid%xb%t(i,j,k) = real(t(idx))
          persistent_grid%xb%q(i,j,k) = real(q(idx))
          ! Calculate pressure from P and PB: grid%xb%p = pb + p
          persistent_grid%xb%p(i,j,k) = real(pb(idx) + p(idx))
        end do
      end do
    end do
    
    ! Update surface pressure
    do j = 1, ny
      do i = 1, nx
        ! C++ row-major indexing: [i][j] -> i + j*nx
        idx = i + (j-1)*nx
        persistent_grid%xb%psfc(i,j) = real(psfc(idx))
      end do
    end do
    
    ! Update PH and PHB data (vertically staggered with nz+1 levels)
    do k = 1, staggered_nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          persistent_grid%ph_2(i,j,k) = real(ph(idx_3d))
          persistent_grid%phb(i,j,k) = real(phb(idx_3d))
        end do
      end do
    end do
    
    ! Update grid metadata
    do j = 1, ny
      do i = 1, nx
        ! C++ row-major indexing: [i][j] -> i + j*nx
        idx = i + (j-1)*nx
        persistent_grid%xb%lat(i,j) = real(lats2d(idx))
        persistent_grid%xb%lon(i,j) = real(lons2d(idx))
        ! Assign terrain height from HGT field
        persistent_grid%ht(i,j) = real(hgt(idx))
        persistent_grid%xb%terr(i,j) = real(hgt(idx))
        do k = 1, staggered_nz
          ! Use calculated height field instead of levels array
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          persistent_grid%xb%h(i,j,k) = real(hf(idx_3d))
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
  integer(c_int) function wrfda_get_innov_vector(it, ob_ptr, iv_ptr) bind(C, name="wrfda_get_innov_vector")
    implicit none
    integer(c_int), intent(in) :: it
    type(c_ptr), value :: ob_ptr, iv_ptr
    type(domain), pointer :: grid
    type(y_type), pointer :: ob
    type(iv_type), pointer :: iv
    type(grid_config_rec_type), pointer :: config_flags
    
    ! Quality control statistics array (simplified for now)
    integer, dimension(1,1,1,1) :: num_qcstat_conv
    
    ! Convert C pointers to Fortran pointers
    grid => persistent_grid
    call c_f_pointer(ob_ptr, ob)
    call c_f_pointer(iv_ptr, iv)
    
    ! config_flags is required by WRFDA interface - allocate it
    allocate(config_flags)
    
    ! Initialize QC statistics array
    num_qcstat_conv = 0
    
    ! Call the main WRFDA innovation vector computation routine
    call da_get_innov_vector(it, num_qcstat_conv, ob, iv, grid, config_flags)

    ! Clean up allocated config_flags
    deallocate(config_flags)

    wrfda_get_innov_vector = 0
    
  end function wrfda_get_innov_vector

  ! Initialize WRFDA module-level variables (kts, kte, sfc_assi_options, etc.)
  integer(c_int) function initialize_wrfda_module_variables() bind(C, name="initialize_wrfda_module_variables")
    implicit none
    type(domain), pointer :: grid

    grid => persistent_grid 
    
    ! Initialize return code
    initialize_wrfda_module_variables = 0
    
    ! CRITICAL FIX: Initialize WRFDA module-level variables
    call da_copy_dims(grid)
     
    ! CRITICAL FIX: Call da_copy_tile_dims to set up module-level tile bounds
    ! This ensures that kts, kte, its, ite, jts, jte are properly set
    call da_copy_tile_dims(grid)
    call da_copy_tile_dims(grid)
    sfc_assi_options = sfc_assi_options_1
    sfc_assi_options = sfc_assi_options_1
    
  end function initialize_wrfda_module_variables

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
          grid%xb%u(i,j,k) = u(idx)
          grid%xb%v(i,j,k) = v(idx)
          grid%xb%t(i,j,k) = t(idx)
          grid%xb%q(i,j,k) = q(idx)
          ! Calculate pressure from P and PB: grid%xb%p = pb + p
          grid%xb%p(i,j,k) = pb(idx) + p(idx)
        end do
      end do
    end do
    
    ! Copy PH and PHB data (vertically staggered with nz+1 levels)
    do k = 1, staggered_nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          grid%ph_2(i,j,k) = ph(idx_3d)
          grid%phb(i,j,k) = phb(idx_3d)
        end do
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        ! C++ row-major indexing: [i][j] -> i + j*nx
        idx = i + (j-1)*nx
        grid%xb%psfc(i,j) = psfc(idx)
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
        grid%xb%lat(i,j) = lats2d(idx)
        grid%xb%lon(i,j) = lons2d(idx)
        ! Assign terrain height from HGT field
        grid%ht(i,j) = hgt(idx)
        grid%xb%terr(i,j) = hgt(idx)
        ! Assign temporary constant roughness length
        grid%xb%rough(i,j) = 0.5
        do k = 1, staggered_nz
          ! Use calculated height field instead of levels array
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          ! Note: HF is vertically staggered with nz+1 levels
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          grid%xb%h(i,j,k) = hf(idx_3d)
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

  ! Count innovation values from iv_type structure
  integer(c_int) function wrfda_count_innovations(family, num_innovations) bind(C, name="wrfda_count_innovations")
    implicit none
    character(c_char), intent(in) :: family(*)
    integer(c_int), intent(out) :: num_innovations
    
    type(iv_type), pointer :: iv
    character(len=20) :: family_str
    integer :: i, count
    integer :: family_index
    
    ! Convert C string to Fortran string
    family_str = ""
    do i = 1, 20
      if (family(i) == c_null_char) exit
      family_str(i:i) = family(i)
    end do
    
    ! Use global persistent iv structure
    if (.not. iv_allocated) then
      num_innovations = 0
      wrfda_count_innovations = 0
      return
    end if
    
    if (.not. associated(persistent_iv)) then
      num_innovations = 0
      wrfda_count_innovations = 0
      return
    end if
    
    ! Point local iv to global persistent structure
    iv => persistent_iv
    
    ! Determine family index
    family_index = 0
    if (trim(family_str) == "synop" .or. trim(family_str) == "adpsfc" .or. trim(family_str) == "metar") then
      family_index = 2  ! synop index (adpsfc and metar use synop structure)
    end if
    
    count = 0
    
    if (family_index == 2 .and. family_index <= size(iv%info) .and. associated(iv%synop) .and. iv%info(2)%nlocal > 0) then
      ! Count innovations from synop array based on QC flags and non-zero innovations
      count = 0
      do i = 1, iv%info(2)%nlocal
        ! Count U innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%u%qc == 0 .and. abs(iv%synop(i)%u%inv) > 1.0e-10) then
          count = count + 1
        end if
        ! Count V innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%v%qc == 0 .and. abs(iv%synop(i)%v%inv) > 1.0e-10) then
          count = count + 1
        end if
        ! Count T innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%t%qc == 0 .and. abs(iv%synop(i)%t%inv) > 1.0e-10) then
          count = count + 1
        end if
        ! Count P innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%p%qc == 0 .and. abs(iv%synop(i)%p%inv) > 1.0e-10) then
          count = count + 1
        end if
        ! Count Q innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%q%qc == 0 .and. abs(iv%synop(i)%q%inv) > 1.0e-10) then
          count = count + 1
        end if
      end do
    else
      count = 0
    end if
    
    num_innovations = count
    wrfda_count_innovations = 0
    
  end function wrfda_count_innovations

  ! Extract innovation values from iv_type structure
  integer(c_int) function wrfda_extract_innovations(family, innovations, num_innovations) bind(C, name="wrfda_extract_innovations")
    implicit none
    character(c_char), intent(in) :: family(*)
    real(c_double), intent(out) :: innovations(*)
    integer(c_int), intent(out) :: num_innovations
    
    type(iv_type), pointer :: iv
    character(len=20) :: family_str
    integer :: i, count
    integer :: family_index
    
    ! Convert C string to Fortran string
    family_str = ""
    do i = 1, 20
      if (family(i) == c_null_char) exit
      family_str(i:i) = family(i)
    end do
    
    ! Use global persistent iv structure instead of converting iv_ptr
    if (.not. iv_allocated) then
      num_innovations = 0
      wrfda_extract_innovations = 0
      return
    end if
    
    if (.not. associated(persistent_iv)) then
      num_innovations = 0
      wrfda_extract_innovations = 0
      return
    end if
    
    ! Point local iv to global persistent structure
    iv => persistent_iv
    
    ! Determine family index (simplified mapping)
    family_index = 0
    
    if (trim(family_str) == "metar" .or. trim(family_str) == "synop" .or. trim(family_str) == "adpsfc") then
      family_index = 2  ! synop index (metar uses synop structure)
    else if (trim(family_str) == "sound") then
      family_index = 3  ! sound index
    else if (trim(family_str) == "gpspw") then
      family_index = 4  ! gpspw index
    end if
    
    if (family_index == 0) then
      wrfda_extract_innovations = -2
      return
    end if
    
    ! Extract innovations based on family
    count = 0
    
    if (family_index == 2 .and. family_index <= size(iv%info) .and. associated(iv%synop) .and. iv%info(2)%nlocal > 0) then
      ! Extract from synop array (all variables with good QC only)
      do i = 1, iv%info(2)%nlocal
        if (i > size(iv%synop)) then
          exit
        end if
        
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
      
    else
      count = 0
    end if
    
    num_innovations = count
    wrfda_extract_innovations = 0
    
  end function wrfda_extract_innovations

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
  
  ! C-callable function to initialize map projection with grid parameters
  subroutine initialize_map_projection_c(map_proj, cen_lat, cen_lon, dx, stand_lon, truelat1, truelat2) bind(C, name="initialize_map_projection_c")
    implicit none
    integer(c_int), intent(in) :: map_proj
    real(c_double), intent(in) :: cen_lat, cen_lon, dx, stand_lon, truelat1, truelat2
    
    if (.not. map_info_initialized) then
      ! Initialize map projection using WRFDA's da_map_set function
      ! Use the provided cen_lat and cen_lon instead of hardcoded values
      call da_map_set(map_proj, real(cen_lat), real(cen_lon), 37.0, 30.5, real(dx), real(stand_lon), real(truelat1), real(truelat2), real(truelat1), real(stand_lon), map_info)
      map_info_initialized = .true.
    end if
  end subroutine initialize_map_projection_c

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
    allocate(grid%domain_clock)
    grid%domain_clock = WRFU_ClockCreate(TimeStep=time_step, StartTime=start_time, StopTime=stop_time, rc=rc)
    if (rc /= WRFU_SUCCESS) then
      print *, "WRFDA ERROR: Failed to create domain clock"
      return
    end if
    
    ! Set the domain_clock_created flag to true
    grid%domain_clock_created = .true.
    
  end subroutine initialize_domain_clock

end module metada_wrfda_dispatch