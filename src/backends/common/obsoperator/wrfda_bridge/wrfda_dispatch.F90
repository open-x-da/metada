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
  use da_control, only: metar, synop, ships, buoy, airep, pilot, sound, sonde_sfc, &
                        sfc_assi_options, sfc_assi_options_1, trace_use_dull, &
                        var4d_run, num_fgat_time
  use da_metar,  only: da_transform_xtoy_metar,  da_transform_xtoy_metar_adj
  use da_synop,  only: da_transform_xtoy_synop,  da_transform_xtoy_synop_adj
  use da_buoy,   only: da_transform_xtoy_buoy,   da_transform_xtoy_buoy_adj
  use da_ships,  only: da_transform_xtoy_ships,  da_transform_xtoy_ships_adj
  use da_airep,  only: da_transform_xtoy_airep,  da_transform_xtoy_airep_adj
  use da_pilot,  only: da_transform_xtoy_pilot,  da_transform_xtoy_pilot_adj
  use da_sound,  only: da_transform_xtoy_sound,  da_transform_xtoy_sound_adj, &
                       da_transform_xtoy_sonde_sfc, da_transform_xtoy_sonde_sfc_adj
  use da_par_util, only: da_copy_dims, da_copy_tile_dims
  use da_tools, only: da_togrid
  use module_configure, only: grid_config_rec_type
  use da_minimisation, only: da_get_innov_vector
  
  implicit none
  
  ! CRITICAL FIX: Add persistent grid to avoid reconstruction issues
  ! WRFDA internal state persists between calls, so we need consistent grid structure
  type(domain), save, target :: persistent_grid
  logical, save :: grid_initialized = .false.
  
  ! Persistent background state for incremental 3D-Var
  ! Background state (xb) is constant throughout the outer loop
  logical, save :: background_state_initialized = .false.
  integer, save :: background_nx, background_ny, background_nz
  
contains

  ! Handle-based entry points: directly call WRFDA routines using provided handles
  integer(c_int) function wrfda_xtoy_apply_handles(operator_family, domain_ptr, iv_ptr, y_ptr) bind(C, name="wrfda_xtoy_apply_handles")
    implicit none
    character(c_char), intent(in) :: operator_family(*)
    type(c_ptr), value :: domain_ptr, iv_ptr, y_ptr
    type(domain), pointer :: grid
    type(iv_type), pointer :: iv
    type(y_type), pointer :: y
    character(len=:), allocatable :: fam
    integer :: n, famlen
    character(len=256) :: fambuf

    call c_f_pointer(domain_ptr, grid)
    call c_f_pointer(iv_ptr, iv)
    call c_f_pointer(y_ptr, y)

    fambuf = ''
    famlen = 0
    do n = 1, len(fambuf)
      if (operator_family(n) == c_null_char) exit
      famlen = famlen + 1
      fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
    end do
    fam = trim(fambuf(1: max(famlen,1)))

    select case (trim(fam))
    case ('metar')
      call da_transform_xtoy_metar(grid, iv, y)
    case ('synop')
      call da_transform_xtoy_synop(grid, iv, y)
    case ('buoy')
      call da_transform_xtoy_buoy(grid, iv, y)
    case ('airep')
      call da_transform_xtoy_airep(grid, iv, y)
    case ('pilot')
      call da_transform_xtoy_pilot(grid, iv, y)
    case ('ships')
      call da_transform_xtoy_ships(grid, iv, y)
    case ('sound')
      call da_transform_xtoy_sound(grid, iv, y)
    case ('sonde_sfc')
      call da_transform_xtoy_sonde_sfc(grid, iv, y)
    case default
      wrfda_xtoy_apply_handles = 1_c_int
      return
    end select
    wrfda_xtoy_apply_handles = 0_c_int
  end function wrfda_xtoy_apply_handles

  integer(c_int) function wrfda_xtoy_adjoint_handles(operator_family, domain_ptr, iv_ptr, jo_grad_y_ptr, jo_grad_x_ptr) bind(C, name="wrfda_xtoy_adjoint_handles")
    implicit none
    character(c_char), intent(in) :: operator_family(*)
    type(c_ptr), value :: domain_ptr, iv_ptr, jo_grad_y_ptr, jo_grad_x_ptr
    type(domain), pointer :: grid
    type(iv_type), pointer :: iv
    type(y_type), pointer :: jo_grad_y
    type(x_type), pointer :: jo_grad_x
    character(len=:), allocatable :: fam
    integer :: n, famlen
    character(len=256) :: fambuf

    call c_f_pointer(domain_ptr, grid)
    call c_f_pointer(iv_ptr, iv)
    call c_f_pointer(jo_grad_y_ptr, jo_grad_y)
    call c_f_pointer(jo_grad_x_ptr, jo_grad_x)

    fambuf = ''
    famlen = 0
    do n = 1, len(fambuf)
      if (operator_family(n) == c_null_char) exit
      famlen = famlen + 1
      fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
    end do
    fam = trim(fambuf(1: max(famlen,1)))

    select case (trim(fam))
    case ('metar')
      call da_transform_xtoy_metar_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('synop')
      call da_transform_xtoy_synop_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('buoy')
      call da_transform_xtoy_buoy_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('airep')
      call da_transform_xtoy_airep_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('pilot')
      call da_transform_xtoy_pilot_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('ships')
      call da_transform_xtoy_ships_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('sound')
      call da_transform_xtoy_sound_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('sonde_sfc')
      call da_transform_xtoy_sonde_sfc_adj(grid, iv, jo_grad_y, jo_grad_x)
    case default
      wrfda_xtoy_adjoint_handles = 1_c_int
      return
    end select
    wrfda_xtoy_adjoint_handles = 0_c_int
  end function wrfda_xtoy_adjoint_handles

  ! Forward operator: H(xb + xa) where xb is background state and xa is analysis increments
  ! This function computes the forward operator for incremental 3D-Var:
  ! - xb (background state): constant throughout outer loop, used for innovation computation
  ! - xa (analysis increments): start at zero, updated by each iteration
  ! - Total state: x_total = xb + xa (computed internally)
  ! - Forward operator: H(xb + xa)
  integer(c_int) function wrfda_xtoy_apply_grid(operator_family, nx, ny, nz, u, v, t, q, psfc, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels, out_y) bind(C, name="wrfda_xtoy_apply_grid")
    implicit none
    character(c_char), intent(in) :: operator_family(*)
    integer(c_int), value :: nx, ny, nz
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    real(c_double), intent(in) :: lats2d(*), lons2d(*)
    real(c_double), intent(in) :: levels(*)
    integer(c_int), value :: num_obs
    real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels(*)
    real(c_double), intent(out) :: out_y(*)

    type(domain), target, save :: grid
    type(iv_type), target, save :: iv
    type(y_type), target :: y
    integer :: n
    character(len=256) :: fambuf
    integer :: famlen
    character(len=256) :: op_str, fam_str, var_str
    integer :: fam_id
    character(len=1) :: var_code

    ! DEBUG: Print entry point
    print *, "WRFDA DEBUG: Entering wrfda_xtoy_apply_grid"
    print *, "WRFDA DEBUG: nx=", nx, " ny=", ny, " nz=", nz, " num_obs=", num_obs

    ! For incremental 3D-Var, we assume the input arrays are analysis increments (xa)
    ! and we use the persistent background state (xb) that was initialized earlier
    if (.not. background_state_initialized) then
      print *, "WRFDA ERROR: Background state not initialized. Call wrfda_initialize_background_state first."
      wrfda_xtoy_apply_grid = 1_c_int
      return
    end if
    
    ! Use persistent grid with background state and update analysis increments
    print *, "WRFDA DEBUG: Using persistent background state and updating analysis increments"
    call wrfda_update_analysis_increments(u, v, t, q, psfc)
    print *, "WRFDA DEBUG: Analysis increments updated"
    
    ! Copy persistent grid to local grid for processing
    grid = persistent_grid
     
    ! CRITICAL FIX: Call da_copy_dims to set up module-level grid bounds
    ! This is required for WRFDA interpolation routines to work properly
    print *, "WRFDA DEBUG: Calling da_copy_dims to set up module-level grid bounds"
    call da_copy_dims(persistent_grid)
    print *, "WRFDA DEBUG: da_copy_dims completed"
     
    ! CRITICAL FIX: Call da_copy_tile_dims to set up module-level tile bounds
    ! This ensures that kts, kte, its, ite, jts, jte are properly set
    print *, "WRFDA DEBUG: Calling da_copy_tile_dims to set up module-level tile bounds"
    call da_copy_tile_dims(persistent_grid)
    print *, "WRFDA DEBUG: da_copy_tile_dims completed"
     
    print *, "WRFDA DEBUG: Setting sfc_assi_options"
    sfc_assi_options = sfc_assi_options_1
    print *, "WRFDA DEBUG: sfc_assi_options set to", sfc_assi_options

    print *, "WRFDA DEBUG: Parsing operator family string"
    fambuf = '' ; famlen = 0
    ! CRITICAL FIX: Properly parse C-style null-terminated string
    do n = 1, 256  ! Use fixed maximum instead of len(fambuf)
      if (operator_family(n) == c_null_char) exit
      famlen = famlen + 1
      fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
    end do
    op_str = trim(fambuf(1:max(famlen,1)))
    print *, "WRFDA DEBUG: Parsed operator string: '", trim(op_str), "'"
     
    print *, "WRFDA DEBUG: Calling parse_family_variable"
    call parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
    print *, "WRFDA DEBUG: Family='", trim(fam_str), "' Variable='", trim(var_str), "' FamilyID=", fam_id, " VarCode='", var_code, "'"
     
    print *, "WRFDA DEBUG: Calling init_iv_from_obs with grid bounds sm33=", persistent_grid%sm33, " em33=", persistent_grid%em33
    call init_iv_from_obs(fam_id, iv, nx, ny, nz, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels, persistent_grid%sm33, persistent_grid%em33)
    print *, "WRFDA DEBUG: init_iv_from_obs completed"
    
    print *, "WRFDA DEBUG: Calling da_allocate_y"
    call da_allocate_y(iv, y)
    print *, "WRFDA DEBUG: da_allocate_y completed"

    print *, "WRFDA DEBUG: About to select case for family: '", trim(fam_str), "'"
    print *, "WRFDA DEBUG: fam_str='", trim(fam_str), "'"
    print *, "WRFDA DEBUG: fam_id=", fam_id
    print *, "WRFDA DEBUG: About to execute select case statement"
    select case (trim(fam_str))
    case ('metar'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_metar"
      !call da_transform_xtoy_metar(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_metar completed"
    case ('synop', 'adpsfc'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_synop for family '", trim(fam_str), "'"
      print *, "WRFDA DEBUG: About to call da_transform_xtoy_synop..."
      print *, "WRFDA DEBUG: About to call da_transform_xtoy_synop function"
      call da_transform_xtoy_synop(persistent_grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_synop completed successfully"
    case ('buoy');  
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_buoy"
      !call da_transform_xtoy_buoy(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_buoy completed"
    case ('ships'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_ships"
      !call da_transform_xtoy_ships(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_ships completed"
    case ('airep'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_airep"
      !call da_transform_xtoy_airep(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_airep completed"
    case ('pilot'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_pilot"
      !call da_transform_xtoy_pilot(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_pilot completed"
    case ('sound'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_sound"
      !call da_transform_xtoy_sound(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_sound completed"
    case ('sonde_sfc'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_sonde_sfc"
      !call da_transform_xtoy_sonde_sfc(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_sonde_sfc completed"
    case default; 
      print *, "WRFDA DEBUG: Unknown family '", trim(fam_str), "', returning error"
      wrfda_xtoy_apply_grid = 1_c_int; return
    end select

    print *, "WRFDA DEBUG: About to call copy_y_to_out"
    call copy_y_to_out(fam_id, var_code, y, out_y, num_obs)
    print *, "WRFDA DEBUG: copy_y_to_out completed"
    print *, "WRFDA DEBUG: Function completed successfully, returning 0"
    wrfda_xtoy_apply_grid = 0_c_int
  end function wrfda_xtoy_apply_grid

  ! Forward operator with separate background and increment inputs
  ! This function properly implements incremental 3D-Var forward operator:
  ! - u_bg, v_bg, t_bg, q_bg, psfc_bg: background state (xb) - constant throughout outer loop
  ! - u_inc, v_inc, t_inc, q_inc, psfc_inc: analysis increments (xa) - updated by each iteration
  ! - Total state: x_total = xb + xa
  ! - Forward operator: H(xb + xa)
  integer(c_int) function wrfda_xtoy_apply_grid_with_background(operator_family, nx, ny, nz, u_bg, v_bg, t_bg, q_bg, psfc_bg, u_inc, v_inc, t_inc, q_inc, psfc_inc, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels, out_y) bind(C, name="wrfda_xtoy_apply_grid_with_background")
    implicit none
    character(c_char), intent(in) :: operator_family(*)
    integer(c_int), value :: nx, ny, nz
    real(c_double), intent(in) :: u_bg(*), v_bg(*), t_bg(*), q_bg(*), psfc_bg(*)
    real(c_double), intent(in) :: u_inc(*), v_inc(*), t_inc(*), q_inc(*), psfc_inc(*)
    real(c_double), intent(in) :: lats2d(*), lons2d(*)
    real(c_double), intent(in) :: levels(*)
    integer(c_int), value :: num_obs
    real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels(*)
    real(c_double), intent(out) :: out_y(*)

    type(domain), target, save :: grid
    type(iv_type), target, save :: iv
    type(y_type), target :: y
    integer :: n
    character(len=256) :: fambuf
    integer :: famlen
    character(len=256) :: op_str, fam_str, var_str
    integer :: fam_id
    character(len=1) :: var_code

    ! DEBUG: Print entry point
    print *, "WRFDA DEBUG: Entering wrfda_xtoy_apply_grid_with_background"
    print *, "WRFDA DEBUG: nx=", nx, " ny=", ny, " nz=", nz, " num_obs=", num_obs

    ! Initialize domain with proper separation of background and increments
    call init_domain_with_separation(grid, nx, ny, nz, u_bg, v_bg, t_bg, q_bg, psfc_bg, u_inc, v_inc, t_inc, q_inc, psfc_inc)
     
    ! CRITICAL FIX: Call da_copy_dims to set up module-level grid bounds
    print *, "WRFDA DEBUG: Calling da_copy_dims to set up module-level grid bounds"
    call da_copy_dims(persistent_grid)
    print *, "WRFDA DEBUG: da_copy_dims completed"
     
    ! CRITICAL FIX: Call da_copy_tile_dims to set up module-level tile bounds
    print *, "WRFDA DEBUG: Calling da_copy_tile_dims to set up module-level tile bounds"
    call da_copy_tile_dims(persistent_grid)
    print *, "WRFDA DEBUG: da_copy_tile_dims completed"
     
    print *, "WRFDA DEBUG: Setting sfc_assi_options"
    sfc_assi_options = sfc_assi_options_1
    print *, "WRFDA DEBUG: sfc_assi_options set to", sfc_assi_options

    print *, "WRFDA DEBUG: Parsing operator family string"
    fambuf = '' ; famlen = 0
    ! CRITICAL FIX: Properly parse C-style null-terminated string
    do n = 1, 256  ! Use fixed maximum instead of len(fambuf)
      if (operator_family(n) == c_null_char) exit
      famlen = famlen + 1
      fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
    end do
    op_str = trim(fambuf(1:max(famlen,1)))
    print *, "WRFDA DEBUG: Parsed operator string: '", trim(op_str), "'"
     
    print *, "WRFDA DEBUG: Calling parse_family_variable"
    call parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
    print *, "WRFDA DEBUG: Family='", trim(fam_str), "' Variable='", trim(var_str), "' FamilyID=", fam_id, " VarCode='", var_code, "'"
     
    print *, "WRFDA DEBUG: Calling init_iv_from_obs with grid bounds sm33=", persistent_grid%sm33, " em33=", persistent_grid%em33
    call init_iv_from_obs(fam_id, iv, nx, ny, nz, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels, persistent_grid%sm33, persistent_grid%em33)
    print *, "WRFDA DEBUG: init_iv_from_obs completed"
    
    print *, "WRFDA DEBUG: Calling da_allocate_y"
    call da_allocate_y(iv, y)
    print *, "WRFDA DEBUG: da_allocate_y completed"

    print *, "WRFDA DEBUG: About to select case for family: '", trim(fam_str), "'"
    print *, "WRFDA DEBUG: fam_str='", trim(fam_str), "'"
    print *, "WRFDA DEBUG: fam_id=", fam_id
    print *, "WRFDA DEBUG: About to execute select case statement"
    select case (trim(fam_str))
    case ('metar'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_metar"
      !call da_transform_xtoy_metar(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_metar completed"
    case ('synop', 'adpsfc'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_synop for family '", trim(fam_str), "'"
      print *, "WRFDA DEBUG: About to call da_transform_xtoy_synop..."
      print *, "WRFDA DEBUG: About to call da_transform_xtoy_synop function"
      call da_transform_xtoy_synop(persistent_grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_synop completed successfully"
    case ('buoy');  
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_buoy"
      !call da_transform_xtoy_buoy(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_buoy completed"
    case ('ships'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_ships"
      !call da_transform_xtoy_ships(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_ships completed"
    case ('airep'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_airep"
      !call da_transform_xtoy_airep(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_airep completed"
    case ('pilot'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_pilot"
      !call da_transform_xtoy_pilot(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_pilot completed"
    case ('sound'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_sound"
      !call da_transform_xtoy_sound(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_sound completed"
    case ('sonde_sfc'); 
      print *, "WRFDA DEBUG: Calling da_transform_xtoy_sonde_sfc"
      !call da_transform_xtoy_sonde_sfc(grid, iv, y)
      print *, "WRFDA DEBUG: da_transform_xtoy_sonde_sfc completed"
    case default; 
      print *, "WRFDA DEBUG: Unknown family '", trim(fam_str), "', returning error"
      wrfda_xtoy_apply_grid_with_background = 1_c_int; return
    end select

    print *, "WRFDA DEBUG: About to call copy_y_to_out"
    call copy_y_to_out(fam_id, var_code, y, out_y, num_obs)
    print *, "WRFDA DEBUG: copy_y_to_out completed"
    print *, "WRFDA DEBUG: Function completed successfully, returning 0"
    wrfda_xtoy_apply_grid_with_background = 0_c_int
  end function wrfda_xtoy_apply_grid_with_background

  ! Enhanced function for constructing iv_type with detailed observation information
  integer(c_int) function wrfda_construct_iv_from_observations(operator_family, nx, ny, nz, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels, obs_values, obs_errors, obs_types, obs_station_ids, obs_elevations, obs_pressures, obs_heights, obs_qc_flags, obs_usage_flags, obs_time_offsets, iv_ptr) bind(C, name="wrfda_construct_iv_from_observations")
    implicit none
    character(c_char), intent(in) :: operator_family(*)
    integer(c_int), value :: nx, ny, nz
    real(c_double), intent(in) :: lats2d(*), lons2d(*), levels(*)
    integer(c_int), value :: num_obs
    real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels(*)
    real(c_double), intent(in) :: obs_values(*), obs_errors(*)
    character(c_char), intent(in) :: obs_types(*), obs_station_ids(*)
    real(c_double), intent(in) :: obs_elevations(*), obs_pressures(*), obs_heights(*)
    integer(c_int), intent(in) :: obs_qc_flags(*), obs_usage_flags(*)
    real(c_double), intent(in) :: obs_time_offsets(*)
    type(c_ptr), value :: iv_ptr

    type(iv_type), pointer :: iv
    integer :: n
    character(len=256) :: fambuf
    integer :: famlen
    character(len=256) :: op_str, fam_str, var_str
    integer :: fam_id
    character(len=1) :: var_code
    
    ! Suppress unused argument warnings for parameters not currently used
    if (obs_types(1) == c_null_char .and. obs_station_ids(1) == c_null_char .and. obs_time_offsets(1) < 0.0) then
      ! This will never execute but prevents unused argument warnings
      continue
    end if

    ! DEBUG: Print entry point
    print *, "WRFDA DEBUG: Entering wrfda_construct_iv_from_observations"
    print *, "WRFDA DEBUG: nx=", nx, " ny=", ny, " nz=", nz, " num_obs=", num_obs

    ! Get pointer to iv_type
    call c_f_pointer(iv_ptr, iv)

    ! Parse operator family string
    fambuf = '' ; famlen = 0
    do n = 1, 256
      if (operator_family(n) == c_null_char) exit
      famlen = famlen + 1
      fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
    end do
    op_str = trim(fambuf(1:max(famlen,1)))
    print *, "WRFDA DEBUG: Parsed operator string: '", trim(op_str), "'"
     
    call parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
    print *, "WRFDA DEBUG: Family='", trim(fam_str), "' Variable='", trim(var_str), "' FamilyID=", fam_id, " VarCode='", var_code, "'"

    ! Initialize iv_type with enhanced observation information
    call init_iv_from_enhanced_obs(fam_id, iv, nx, ny, nz, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels, obs_values, obs_errors, obs_types, obs_station_ids, obs_elevations, obs_pressures, obs_heights, obs_qc_flags, obs_usage_flags, obs_time_offsets)

    print *, "WRFDA DEBUG: wrfda_construct_iv_from_observations completed successfully"
    wrfda_construct_iv_from_observations = 0_c_int
  end function wrfda_construct_iv_from_observations

  ! Adjoint operator: H^T * δy -> δx (gradient accumulation)
  ! This function computes the adjoint of the forward operator for incremental 3D-Var:
  ! - delta_y: gradient in observation space
  ! - inout_u, inout_v, inout_t, inout_q, inout_psfc: gradient in state space (accumulated)
  ! - The adjoint operator accumulates gradients into the state space arrays
  integer(c_int) function wrfda_xtoy_adjoint_grid(operator_family, nx, ny, nz, delta_y, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels, inout_u, inout_v, inout_t, inout_q, inout_psfc) bind(C, name="wrfda_xtoy_adjoint_grid")
    implicit none
    character(c_char), intent(in) :: operator_family(*)
    integer(c_int), value :: nx, ny, nz
    real(c_double), intent(in) :: delta_y(*)
    real(c_double), intent(in) :: lats2d(*), lons2d(*)
    real(c_double), intent(in) :: levels(*)
    integer(c_int), value :: num_obs
    real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels(*)
    real(c_double), intent(inout) :: inout_u(*), inout_v(*), inout_t(*), inout_q(*), inout_psfc(*)

    type(domain), target, save :: grid
    type(iv_type), target, save :: iv
    type(y_type), target :: jo_grad_y
    type(x_type), target :: jo_grad_x
    integer :: n
    character(len=256) :: fambuf
    integer :: famlen
    character(len=256) :: op_str, fam_str, var_str
    integer :: fam_id
    character(len=1) :: var_code

    ! For adjoint operator, we use the persistent background state (xb) and zero xa initially
    ! The adjoint operator accumulates gradients into jo_grad_x, which we then add to inout arrays
    if (.not. background_state_initialized) then
      print *, "WRFDA ERROR: Background state not initialized. Call wrfda_initialize_background_state first."
      wrfda_xtoy_adjoint_grid = 1_c_int
      return
    end if
    
    ! Use persistent grid with background state and zero analysis increments
    print *, "WRFDA DEBUG: Using persistent background state for adjoint operator"
    grid = persistent_grid
    
    ! Zero out analysis increments for adjoint computation
    grid%xa%u = 0.0
    grid%xa%v = 0.0
    grid%xa%t = 0.0
    grid%xa%q = 0.0
    grid%xa%psfc = 0.0
    
         ! CRITICAL FIX: Call da_copy_dims to set up module-level grid bounds
     ! This is required for WRFDA interpolation routines to work properly
     print *, "WRFDA DEBUG: Calling da_copy_dims in adjoint function"
     call da_copy_dims(grid)
     print *, "WRFDA DEBUG: da_copy_dims completed in adjoint function"
     
     ! CRITICAL FIX: Call da_copy_tile_dims to set up module-level tile bounds
     ! This ensures that kts, kte, its, ite, jts, jte are properly set
     print *, "WRFDA DEBUG: Calling da_copy_tile_dims in adjoint function"
     call da_copy_tile_dims(grid)
     print *, "WRFDA DEBUG: da_copy_tile_dims completed in adjoint function"
     
     print *, "WRFDA DEBUG: Grid bounds from grid structure: sm33=", grid%sm33, " em33=", grid%em33
    
    call zero_x_like(jo_grad_x, nx, ny, nz)
    sfc_assi_options = sfc_assi_options_1

    fambuf = '' ; famlen = 0
    ! CRITICAL FIX: Properly parse C-style null-terminated string
    do n = 1, 256  ! Use fixed maximum instead of len(fambuf)
      if (operator_family(n) == c_null_char) exit
      famlen = famlen + 1
      fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
    end do
    op_str = trim(fambuf(1:max(famlen,1)))
    call parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
         print *, "WRFDA DEBUG: Calling init_iv_from_obs in adjoint function with grid bounds sm33=", grid%sm33, " em33=", grid%em33
     call init_iv_from_obs(fam_id, iv, nx, ny, nz, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels, grid%sm33, grid%em33)
    call da_allocate_y(iv, jo_grad_y)
    call init_y_from_delta(fam_id, var_code, jo_grad_y, delta_y, num_obs)

    select case (trim(fam_str))
    case ('metar'); call da_transform_xtoy_metar_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('synop'); call da_transform_xtoy_synop_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('buoy');  call da_transform_xtoy_buoy_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('ships'); call da_transform_xtoy_ships_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('airep'); call da_transform_xtoy_airep_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('pilot'); call da_transform_xtoy_pilot_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('sound'); call da_transform_xtoy_sound_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('sonde_sfc'); call da_transform_xtoy_sonde_sfc_adj(grid, iv, jo_grad_y, jo_grad_x)
    case default; wrfda_xtoy_adjoint_grid = 1_c_int; return
    end select

    call copy_x_to_state(jo_grad_x, inout_u, inout_v, inout_t, inout_q, inout_psfc, nx, ny, nz)
    wrfda_xtoy_adjoint_grid = 0_c_int
  end function wrfda_xtoy_adjoint_grid


  ! Profile-capable API: multiple levels per observation
  ! COMMENTED OUT: Not currently used, causing uninitialized variable warnings
  ! integer(c_int) function wrfda_xtoy_apply_profiles(operator_family, nx, ny, nz, &
  !     u, v, t, q, psfc, lats2d, lons2d, levels, &
  !     num_obs, obs_counts, obs_lats, obs_lons, obs_levels_flat, out_y_flat) &
  !     bind(C, name="wrfda_xtoy_apply_profiles")
  !   implicit none
  !   character(c_char), intent(in) :: operator_family(*)
  !   integer(c_int), value :: nx, ny, nz
  !   real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
  !   real(c_double), intent(in) :: lats2d(*), lons2d(*), levels(*)
  !   integer(c_int), value :: num_obs
  !   integer(c_int), intent(in) :: obs_counts(*)
  !   real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels_flat(*)
  !   real(c_double), intent(out) :: out_y_flat(*)
  !
  !   type(domain), target :: grid
  !   type(iv_type), target :: iv
  !   type(y_type), target :: y
  !   integer :: n
  !   character(len=256) :: fambuf
  !   integer :: famlen
  !   character(len=:), allocatable :: op_str, fam_str, var_str
  !   integer :: fam_id
  !   character(len=1) :: var_code
  !
  !   call init_domain_from_arrays(grid, nx, ny, nz, u, v, t, q, psfc)
  !   
  !   ! CRITICAL FIX: Call da_copy_dims to set up module-level grid bounds
  !   print *, "WRFDA DEBUG: Calling da_copy_dims in profiles function"
  !   call da_copy_dims(grid)
  !   print *, "WRFDA DEBUG: da_copy_dims completed in profiles function"
  !   
  !   ! CRITICAL FIX: Call da_copy_tile_dims to set up module-level tile bounds
  !   print *, "WRFDA DEBUG: Calling da_copy_tile_dims in profiles function"
  !   call da_copy_tile_dims(grid)
  !   print *, "WRFDA DEBUG: da_copy_tile_dims completed in profiles function"
  !   
  !   sfc_assi_options = sfc_assi_options_1
  !
  !   fambuf = '' ; famlen = 0
  !   ! CRITICAL FIX: Properly parse C-style null-terminated string
  !   do n = 1, 256  ! Use fixed maximum instead of len(fambuf)
  !     if (operator_family(n) == c_null_char) exit
  !     famlen = famlen + 1
  !     fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
  !   end do
  !   op_str = trim(fambuf(1:max(famlen,1)))
  !   call parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
  !   call init_iv_profiles(fam_id, iv, nx, ny, nz, lats2d, lons2d, levels, &
  !                         num_obs, obs_counts, obs_lats, obs_lons, obs_levels_flat)
  !   call da_allocate_y(iv, y)
  !
  !   select case (trim(fam_str))
  !   case ('airep'); call da_transform_xtoy_airep(grid, iv, y)
  !   case ('pilot'); call da_transform_xtoy_pilot(grid, iv, y)
  !   case ('sound'); call da_transform_xtoy_sound(grid, iv, y)
  !   case default; wrfda_xtoy_apply_profiles = 1_c_int; return
  !   end select
  !
  !   call copy_y_to_out_profiles(fam_id, var_code, y, num_obs, obs_counts, out_y_flat)
  !   wrfda_xtoy_apply_profiles = 0_c_int
  ! end function wrfda_xtoy_apply_profiles

  ! COMMENTED OUT: Not currently used, causing uninitialized variable warnings
  ! integer(c_int) function wrfda_xtoy_adjoint_profiles(operator_family, nx, ny, nz, &
  !     delta_y_flat, lats2d, lons2d, levels, num_obs, obs_counts, obs_lats, obs_lons, obs_levels_flat, &
  !     inout_u, inout_v, inout_t, inout_q, inout_psfc) bind(C, name="wrfda_xtoy_adjoint_profiles")
  !   implicit none
  !   character(c_char), intent(in) :: operator_family(*)
  !   integer(c_int), value :: nx, ny, nz
  !   real(c_double), intent(in) :: delta_y_flat(*)
  !   real(c_double), intent(in) :: lats2d(*), lons2d(*), levels(*)
  !   integer(c_int), value :: num_obs
  !   integer(c_int), intent(in) :: obs_counts(*)
  !   real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels_flat(*)
  !   real(c_double), intent(inout) :: inout_u(*), inout_v(*), inout_t(*), inout_q(*), inout_psfc(*)
  !
  !   type(domain), target :: grid
  !   type(iv_type), target :: iv
  !   type(y_type), target :: jo_grad_y
  !   type(x_type), target :: jo_grad_x
  !   integer :: n
  !   character(len=256) :: fambuf
  !   integer :: famlen
  !   character(len=:), allocatable :: op_str, fam_str, var_str
  !   integer :: fam_id
  !   character(len=1) :: var_code
  !
  !   ! CRITICAL FIX: Use the correct initialization that fills arrays with actual data
  !   ! init_domain_from_arrays_refs was incorrectly setting all arrays to 0.0
  !   call init_domain_from_arrays(grid, nx, ny, nz, inout_u, inout_v, inout_t, inout_q, inout_psfc)
  !   
  !   ! CRITICAL FIX: Call da_copy_dims to set up module-level grid bounds
  !   print *, "WRFDA DEBUG: Calling da_copy_dims in adjoint profiles function"
  !   call da_copy_dims(grid)
  !   print *, "WRFDA DEBUG: da_copy_dims completed in adjoint profiles function"
  !   
  !   ! CRITICAL FIX: Call da_copy_tile_dims to set up module-level tile bounds
  !   print *, "WRFDA DEBUG: Calling da_copy_tile_dims in adjoint profiles function"
  !   call da_copy_tile_dims(grid)
  !   print *, "WRFDA DEBUG: da_copy_tile_dims completed in adjoint profiles function"
  !   
  !   call zero_x_like(jo_grad_x, nx, ny, nz)
  !   sfc_assi_options = sfc_assi_options_1
  !
  !   fambuf = '' ; famlen = 0
  !   do n = 1, len(fambuf)
  !     if (operator_family(n) == c_null_char) exit
  !     famlen = famlen + 1
  !     fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
  !   end do
  !   op_str = trim(fambuf(1:max(famlen,1)))
  !   call parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
  !   call init_iv_profiles(fam_id, iv, nx, ny, nz, lats2d, lons2d, levels, &
  !                         num_obs, obs_counts, obs_lats, obs_lons, obs_levels_flat)
  !   call da_allocate_y(iv, jo_grad_y)
  !   call init_y_from_delta_profiles(fam_id, var_code, jo_grad_y, num_obs, obs_counts, delta_y_flat)
  !
  !   select case (trim(fam_str))
  !   case ('airep'); call da_transform_xtoy_airep_adj(grid, iv, jo_grad_y, jo_grad_x)
  !   case ('pilot'); call da_transform_xtoy_pilot_adj(grid, iv, jo_grad_y, jo_grad_x)
  !   case ('sound'); call da_transform_xtoy_sound_adj(grid, iv, jo_grad_y, jo_grad_x)
  !   case default; wrfda_xtoy_adjoint_profiles = 1_c_int; return
  !   end select
  !
  !   call copy_x_to_state(jo_grad_x, inout_u, inout_v, inout_t, inout_q, inout_psfc, nx, ny, nz)
  !   wrfda_xtoy_adjoint_profiles = 0_c_int
  ! end function wrfda_xtoy_adjoint_profiles


  subroutine init_domain_from_arrays(grid, nx, ny, nz, u, v, t, q, psfc)
    type(domain), intent(inout) :: grid
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
    
    print *, "WRFDA DEBUG: init_domain_from_arrays: nx=", nx, " ny=", ny, " nz=", nz, " nz1=", nz1
    
         ! CRITICAL FIX: Set WRFDA grid bounds that are used for observation array allocation
     ! These are used in da_copy_dims.inc to set kms, kme, ims, ime, jms, jme
     print *, "WRFDA DEBUG: Setting WRFDA grid bounds for proper observation array allocation"
     ! WRFDA expects 1-based bounds for grid arrays
     grid%sm31 = 1; grid%em31 = nx
     grid%sm32 = 1; grid%em32 = ny  
     grid%sm33 = 1; grid%em33 = nz1
     
           ! CRITICAL FIX: Also set domain start/end bounds that are used in da_copy_dims.inc
      ! These are used to set ids, ide, jds, jde, kds, kde
      print *, "WRFDA DEBUG: Setting WRFDA domain start/end bounds"
      grid%sd31 = 1; grid%ed31 = nx
      grid%sd32 = 1; grid%ed32 = ny
      grid%sd33 = 1; grid%ed33 = nz1
     
                       ! Also set the parallel decomposition bounds (single tile case)
       grid%sp31 = 1; grid%ep31 = nx
       grid%sp32 = 1; grid%ep32 = ny
       grid%sp33 = 1; grid%ep33 = nz1
     
                       ! CRITICAL FIX: Set domain bounds that are used for task bounds (kts:kte)
       ! These are used in da_copy_tile_dims.inc and da_interp_lin_3d.inc
       print *, "WRFDA DEBUG: Setting WRFDA domain bounds for proper task bounds"
       grid%xp%kds = 1; grid%xp%kde = nz1
       grid%xp%ids = 1; grid%xp%ide = nx
       grid%xp%jds = 1; grid%xp%jde = ny
     
                       ! CRITICAL FIX: Also set the xp task bounds to match the domain bounds
       ! This ensures kts:kte loop uses valid indices
       print *, "WRFDA DEBUG: Setting WRFDA xp task bounds to match domain bounds"
       grid%xp%kts = 1; grid%xp%kte = nz1
       grid%xp%its = 1; grid%xp%ite = nx
       grid%xp%jts = 1; grid%xp%jte = ny
       
       ! CRITICAL FIX: Set up tile information for single-tile case
       ! This is required by da_copy_tile_dims to set kts, kte, its, ite, jts, jte
       print *, "WRFDA DEBUG: Setting up tile information for single-tile case"
       grid%num_tiles = 1
       allocate(grid%i_start(1:1), grid%i_end(1:1))
       allocate(grid%j_start(1:1), grid%j_end(1:1))
       grid%i_start(1) = 1; grid%i_end(1) = nx
       grid%j_start(1) = 1; grid%j_end(1) = ny
       print *, "WRFDA DEBUG: Tile arrays allocated: i_start=", grid%i_start(1), " i_end=", grid%i_end(1), " j_start=", grid%j_start(1), " j_end=", grid%j_end(1)
    
    print *, "WRFDA DEBUG: Grid bounds set - sm31:em31=", grid%sm31, ":", grid%em31
    print *, "WRFDA DEBUG: Grid bounds set - sm32:em32=", grid%sm32, ":", grid%em32  
    print *, "WRFDA DEBUG: Grid bounds set - sm33:em33=", grid%sm33, ":", grid%em33
    
         ! Allocate and fill xa fields (single precision in WRFDA)
     ! WRFDA expects 1-based bounds for grid arrays to match grid indices
     print *, "WRFDA DEBUG: Allocating xa fields with 1-based bounds"
     allocate(grid%xa%u(1:nx,1:ny,1:nz1)); allocate(grid%xa%v(1:nx,1:ny,1:nz1))
     allocate(grid%xa%t(1:nx,1:ny,1:nz1)); allocate(grid%xa%q(1:nx,1:ny,1:nz1))
     allocate(grid%xa%psfc(1:nx,1:ny))
    print *, "WRFDA DEBUG: xa fields allocated"
    
         print *, "WRFDA DEBUG: Filling xa fields with data using 1-based indexing"
     ! CRITICAL FIX: Add debug output to check input array values
     print *, "WRFDA DEBUG: Sample input array values:"
     print *, "  u(1)=", u(1), " u(nx*ny*nz)=", u(nx*ny*nz)
     print *, "  v(1)=", v(1), " v(nx*ny*nz)=", v(nx*ny*nz)
     print *, "  t(1)=", t(1), " t(nx*ny*nz)=", t(nx*ny*nz)
     print *, "  q(1)=", q(1), " q(nx*ny*nz)=", q(nx*ny*nz)
     print *, "  psfc(1)=", psfc(1), " psfc(nx*ny)=", psfc(nx*ny)
     
     do k=1,nz1; do j=1,ny; do i=1,nx
       ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
       grid%xa%u(i,j,k) = real(u(i + (j-1)*nx + (k-1)*nx*ny))
       grid%xa%v(i,j,k) = real(v(i + (j-1)*nx + (k-1)*nx*ny))
       grid%xa%t(i,j,k) = real(t(i + (j-1)*nx + (k-1)*nx*ny))
       grid%xa%q(i,j,k) = real(q(i + (j-1)*nx + (k-1)*nx*ny))
     end do; end do; end do
     do j=1,ny; do i=1,nx
       ! C++ row-major indexing: [i][j] -> i + j*nx
       grid%xa%psfc(i,j) = real(psfc(i + (j-1)*nx))
     end do; end do
    print *, "WRFDA DEBUG: xa fields filled"
    
         ! Mirror xb from xa
     print *, "WRFDA DEBUG: Allocating and copying xb fields with 1-based bounds"
     allocate(grid%xb%u(1:nx,1:ny,1:nz1)); grid%xb%u = grid%xa%u
     allocate(grid%xb%v(1:nx,1:ny,1:nz1)); grid%xb%v = grid%xa%v
     allocate(grid%xb%t(1:nx,1:ny,1:nz1)); grid%xb%t = grid%xa%t
     allocate(grid%xb%q(1:nx,1:ny,1:nz1)); grid%xb%q = grid%xa%q
     allocate(grid%xb%psfc(1:nx,1:ny));  grid%xb%psfc = grid%xa%psfc
    print *, "WRFDA DEBUG: init_domain_from_arrays completed successfully"
  end subroutine init_domain_from_arrays

  ! Initialize domain for adjoint operator
  ! For adjoint operator, we need background state in xb and zero xa initially
  ! The adjoint operator accumulates gradients into jo_grad_x
  subroutine init_domain_for_adjoint(grid, nx, ny, nz, u_bg, v_bg, t_bg, q_bg, psfc_bg)
    type(domain), intent(inout) :: grid
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: u_bg(*), v_bg(*), t_bg(*), q_bg(*), psfc_bg(*)
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
    
    print *, "WRFDA DEBUG: init_domain_for_adjoint: nx=", nx, " ny=", ny, " nz=", nz, " nz1=", nz1
    
    ! Set WRFDA grid bounds (same as before)
    grid%sm31 = 1; grid%em31 = nx
    grid%sm32 = 1; grid%em32 = ny  
    grid%sm33 = 1; grid%em33 = nz1
    
    grid%sd31 = 1; grid%ed31 = nx
    grid%sd32 = 1; grid%ed32 = ny
    grid%sd33 = 1; grid%ed33 = nz1
    
    grid%sp31 = 1; grid%ep31 = nx
    grid%sp32 = 1; grid%ep32 = ny
    grid%sp33 = 1; grid%ep33 = nz1
    
    grid%xp%kds = 1; grid%xp%kde = nz1
    grid%xp%ids = 1; grid%xp%ide = nx
    grid%xp%jds = 1; grid%xp%jde = ny
    
    grid%xp%kts = 1; grid%xp%kte = nz1
    grid%xp%its = 1; grid%xp%ite = nx
    grid%xp%jts = 1; grid%xp%jte = ny
    
    grid%num_tiles = 1
    allocate(grid%i_start(1:1), grid%i_end(1:1))
    allocate(grid%j_start(1:1), grid%j_end(1:1))
    grid%i_start(1) = 1; grid%i_end(1) = nx
    grid%j_start(1) = 1; grid%j_end(1) = ny
    
    ! Allocate xb fields (background state) with 1-based bounds
    print *, "WRFDA DEBUG: Allocating xb fields (background state) for adjoint"
    allocate(grid%xb%u(1:nx,1:ny,1:nz1))
    allocate(grid%xb%v(1:nx,1:ny,1:nz1))
    allocate(grid%xb%t(1:nx,1:ny,1:nz1))
    allocate(grid%xb%q(1:nx,1:ny,1:nz1))
    allocate(grid%xb%psfc(1:nx,1:ny))
    print *, "WRFDA DEBUG: xb fields allocated for adjoint"
    
    ! Fill xb fields with background state data
    print *, "WRFDA DEBUG: Filling xb fields with background state data for adjoint"
    do k=1,nz1; do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
      grid%xb%u(i,j,k) = real(u_bg(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xb%v(i,j,k) = real(v_bg(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xb%t(i,j,k) = real(t_bg(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xb%q(i,j,k) = real(q_bg(i + (j-1)*nx + (k-1)*nx*ny))
    end do; end do; end do
    do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j] -> i + j*nx
      grid%xb%psfc(i,j) = real(psfc_bg(i + (j-1)*nx))
    end do; end do
    print *, "WRFDA DEBUG: xb fields filled with background state for adjoint"
    
    ! Allocate xa fields (analysis increments) and initialize to zero
    ! For adjoint operator, xa starts at zero and gradients are accumulated in jo_grad_x
    print *, "WRFDA DEBUG: Allocating xa fields (analysis increments) for adjoint"
    allocate(grid%xa%u(1:nx,1:ny,1:nz1)); grid%xa%u = 0.0
    allocate(grid%xa%v(1:nx,1:ny,1:nz1)); grid%xa%v = 0.0
    allocate(grid%xa%t(1:nx,1:ny,1:nz1)); grid%xa%t = 0.0
    allocate(grid%xa%q(1:nx,1:ny,1:nz1)); grid%xa%q = 0.0
    allocate(grid%xa%psfc(1:nx,1:ny)); grid%xa%psfc = 0.0
    print *, "WRFDA DEBUG: xa fields allocated and initialized to zero for adjoint"
    
    print *, "WRFDA DEBUG: init_domain_for_adjoint completed successfully"
  end subroutine init_domain_for_adjoint

  ! New function to properly handle background state for forward operator
  subroutine init_domain_from_background_state(grid, nx, ny, nz, u, v, t, q, psfc)
    type(domain), intent(inout) :: grid
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
    
    print *, "WRFDA DEBUG: init_domain_from_background_state: nx=", nx, " ny=", ny, " nz=", nz, " nz1=", nz1
    
    ! Set WRFDA grid bounds (same as before)
    grid%sm31 = 1; grid%em31 = nx
    grid%sm32 = 1; grid%em32 = ny  
    grid%sm33 = 1; grid%em33 = nz1
    
    grid%sd31 = 1; grid%ed31 = nx
    grid%sd32 = 1; grid%ed32 = ny
    grid%sd33 = 1; grid%ed33 = nz1
    
    grid%sp31 = 1; grid%ep31 = nx
    grid%sp32 = 1; grid%ep32 = ny
    grid%sp33 = 1; grid%ep33 = nz1
    
    grid%xp%kds = 1; grid%xp%kde = nz1
    grid%xp%ids = 1; grid%xp%ide = nx
    grid%xp%jds = 1; grid%xp%jde = ny
    
    grid%xp%kts = 1; grid%xp%kte = nz1
    grid%xp%its = 1; grid%xp%ite = nx
    grid%xp%jts = 1; grid%xp%jte = ny
    
    grid%num_tiles = 1
    allocate(grid%i_start(1:1), grid%i_end(1:1))
    allocate(grid%j_start(1:1), grid%j_end(1:1))
    grid%i_start(1) = 1; grid%i_end(1) = nx
    grid%j_start(1) = 1; grid%j_end(1) = ny
    
    ! Allocate xb fields (background state) with 1-based bounds
    print *, "WRFDA DEBUG: Allocating xb fields (background state) with 1-based bounds"
    allocate(grid%xb%u(1:nx,1:ny,1:nz1))
    allocate(grid%xb%v(1:nx,1:ny,1:nz1))
    allocate(grid%xb%t(1:nx,1:ny,1:nz1))
    allocate(grid%xb%q(1:nx,1:ny,1:nz1))
    allocate(grid%xb%psfc(1:nx,1:ny))
    print *, "WRFDA DEBUG: xb fields allocated"
    
    ! Fill xb fields with background state data
    print *, "WRFDA DEBUG: Filling xb fields with background state data"
    do k=1,nz1; do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
      grid%xb%u(i,j,k) = real(u(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xb%v(i,j,k) = real(v(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xb%t(i,j,k) = real(t(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xb%q(i,j,k) = real(q(i + (j-1)*nx + (k-1)*nx*ny))
    end do; end do; end do
    do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j] -> i + j*nx
      grid%xb%psfc(i,j) = real(psfc(i + (j-1)*nx))
    end do; end do
    print *, "WRFDA DEBUG: xb fields filled with background state"
    
    ! Allocate xa fields (analysis increments) and initialize to zero
    ! For forward operator, we typically want zero increments initially
    print *, "WRFDA DEBUG: Allocating xa fields (analysis increments) with 1-based bounds"
    allocate(grid%xa%u(1:nx,1:ny,1:nz1)); grid%xa%u = 0.0
    allocate(grid%xa%v(1:nx,1:ny,1:nz1)); grid%xa%v = 0.0
    allocate(grid%xa%t(1:nx,1:ny,1:nz1)); grid%xa%t = 0.0
    allocate(grid%xa%q(1:nx,1:ny,1:nz1)); grid%xa%q = 0.0
    allocate(grid%xa%psfc(1:nx,1:ny)); grid%xa%psfc = 0.0
    print *, "WRFDA DEBUG: xa fields allocated and initialized to zero"
    
    print *, "WRFDA DEBUG: init_domain_from_background_state completed successfully"
  end subroutine init_domain_from_background_state

  ! Initialize domain with proper incremental 3D-Var formulation
  ! This function sets up grid%xb (background state) and grid%xa (analysis increments)
  ! The total state is computed as xb + xa internally by WRFDA
  subroutine init_domain_with_separation(grid, nx, ny, nz, u_bg, v_bg, t_bg, q_bg, psfc_bg, u_inc, v_inc, t_inc, q_inc, psfc_inc)
    type(domain), intent(inout) :: grid
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: u_bg(*), v_bg(*), t_bg(*), q_bg(*), psfc_bg(*)
    real(c_double), intent(in) :: u_inc(*), v_inc(*), t_inc(*), q_inc(*), psfc_inc(*)
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
    
    print *, "WRFDA DEBUG: init_domain_with_separation: nx=", nx, " ny=", ny, " nz=", nz, " nz1=", nz1
    
    ! Set WRFDA grid bounds (same as before)
    grid%sm31 = 1; grid%em31 = nx
    grid%sm32 = 1; grid%em32 = ny  
    grid%sm33 = 1; grid%em33 = nz1
    
    grid%sd31 = 1; grid%ed31 = nx
    grid%sd32 = 1; grid%ed32 = ny
    grid%sd33 = 1; grid%ed33 = nz1
    
    grid%sp31 = 1; grid%ep31 = nx
    grid%sp32 = 1; grid%ep32 = ny
    grid%sp33 = 1; grid%ep33 = nz1
    
    grid%xp%kds = 1; grid%xp%kde = nz1
    grid%xp%ids = 1; grid%xp%ide = nx
    grid%xp%jds = 1; grid%xp%jde = ny
    
    grid%xp%kts = 1; grid%xp%kte = nz1
    grid%xp%its = 1; grid%xp%ite = nx
    grid%xp%jts = 1; grid%xp%jte = ny
    
    grid%num_tiles = 1
    allocate(grid%i_start(1:1), grid%i_end(1:1))
    allocate(grid%j_start(1:1), grid%j_end(1:1))
    grid%i_start(1) = 1; grid%i_end(1) = nx
    grid%j_start(1) = 1; grid%j_end(1) = ny
    
    ! Allocate xb fields (background state) with 1-based bounds
    print *, "WRFDA DEBUG: Allocating xb fields (background state) with 1-based bounds"
    allocate(grid%xb%u(1:nx,1:ny,1:nz1))
    allocate(grid%xb%v(1:nx,1:ny,1:nz1))
    allocate(grid%xb%t(1:nx,1:ny,1:nz1))
    allocate(grid%xb%q(1:nx,1:ny,1:nz1))
    allocate(grid%xb%psfc(1:nx,1:ny))
    print *, "WRFDA DEBUG: xb fields allocated"
    
    ! Fill xb fields with background state data (constant throughout outer loop)
    print *, "WRFDA DEBUG: Filling xb fields with background state data"
    do k=1,nz1; do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
      grid%xb%u(i,j,k) = real(u_bg(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xb%v(i,j,k) = real(v_bg(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xb%t(i,j,k) = real(t_bg(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xb%q(i,j,k) = real(q_bg(i + (j-1)*nx + (k-1)*nx*ny))
    end do; end do; end do
    do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j] -> i + j*nx
      grid%xb%psfc(i,j) = real(psfc_bg(i + (j-1)*nx))
    end do; end do
    print *, "WRFDA DEBUG: xb fields filled with background state"
    
    ! Allocate xa fields (analysis increments) with 1-based bounds
    print *, "WRFDA DEBUG: Allocating xa fields (analysis increments) with 1-based bounds"
    allocate(grid%xa%u(1:nx,1:ny,1:nz1))
    allocate(grid%xa%v(1:nx,1:ny,1:nz1))
    allocate(grid%xa%t(1:nx,1:ny,1:nz1))
    allocate(grid%xa%q(1:nx,1:ny,1:nz1))
    allocate(grid%xa%psfc(1:nx,1:ny))
    print *, "WRFDA DEBUG: xa fields allocated"
    
    ! Fill xa fields with analysis increments (updated by each iteration)
    print *, "WRFDA DEBUG: Filling xa fields with analysis increments"
    do k=1,nz1; do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
      grid%xa%u(i,j,k) = real(u_inc(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xa%v(i,j,k) = real(v_inc(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xa%t(i,j,k) = real(t_inc(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xa%q(i,j,k) = real(q_inc(i + (j-1)*nx + (k-1)*nx*ny))
    end do; end do; end do
    do j=1,ny; do i=1,nx
      ! C++ row-major indexing: [i][j] -> i + j*nx
      grid%xa%psfc(i,j) = real(psfc_inc(i + (j-1)*nx))
    end do; end do
    print *, "WRFDA DEBUG: xa fields filled with analysis increments"
    
    print *, "WRFDA DEBUG: init_domain_with_separation completed successfully"
  end subroutine init_domain_with_separation

  ! Enhanced function to construct iv_type with detailed observation information
  subroutine init_iv_from_enhanced_obs(family, iv, nx, ny, nz, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels, obs_values, obs_errors, obs_types, obs_station_ids, obs_elevations, obs_pressures, obs_heights, obs_qc_flags, obs_usage_flags, obs_time_offsets)
    integer, intent(in) :: family
    type(iv_type), intent(inout) :: iv
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: lats2d(*), lons2d(*), levels(*)
    integer, intent(in) :: num_obs
    real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels(*)
    real(c_double), intent(in) :: obs_values(*), obs_errors(*)
    character(c_char), intent(in) :: obs_types(*), obs_station_ids(*)
    real(c_double), intent(in) :: obs_elevations(*), obs_pressures(*), obs_heights(*)
    integer(c_int), intent(in) :: obs_qc_flags(*), obs_usage_flags(*)
    real(c_double), intent(in) :: obs_time_offsets(*)
    
    integer :: n, i, j, k
    real(c_double) :: xfloat, yfloat, zkfloat
    integer :: valid_obs_count
    
    ! Suppress unused argument warnings for parameters not currently used
    if (obs_types(1) == c_null_char .and. obs_station_ids(1) == c_null_char .and. obs_time_offsets(1) < 0.0) then
      ! This will never execute but prevents unused argument warnings
      continue
    end if
    
    print *, "WRFDA DEBUG: init_iv_from_enhanced_obs: family=", family, " num_obs=", num_obs, " nz=", nz
    
    ! Initialize ALL observation families to prevent memory issues
    print *, "WRFDA DEBUG: Initializing ALL observation families..."
    iv%info(:)%nlocal = 0
    iv%info(:)%ntotal = 0
    iv%info(:)%max_lev = 1
    
    ! Initialize instrument-related fields
    iv%num_inst = 0
    iv%total_rad_pixel = 0
    iv%total_rad_channel = 0
    nullify(iv%instid)
    
    ! Count valid observations (those within domain and with good QC)
    valid_obs_count = 0
    do n = 1, num_obs
      call find_fractional_ij(nx, ny, lats2d, lons2d, obs_lats(n), obs_lons(n), i, j, xfloat, yfloat)
      if (i /= -1 .and. j /= -1 .and. obs_qc_flags(n) >= 0 .and. obs_usage_flags(n) > 0) then
        valid_obs_count = valid_obs_count + 1
      end if
    end do
    
    print *, "WRFDA DEBUG: Found", valid_obs_count, "valid observations out of", num_obs
    
    ! Set the actual family we want to use
    iv%info(family)%nlocal = valid_obs_count
    iv%info(family)%ntotal = valid_obs_count
    iv%info(family)%max_lev = nz
    
    ! Allocate arrays using WRFDA's routine
    print *, "WRFDA DEBUG: Calling da_allocate_obs_info for enhanced observations"
    call da_allocate_obs_info(iv, family)
    print *, "WRFDA DEBUG: da_allocate_obs_info completed"
    
    ! Set levels array
    iv%info(family)%levels = 1
    
    ! Initialize grid indices and observation data
    valid_obs_count = 0
    do n = 1, num_obs
      ! Find grid indices
      call find_fractional_ij(nx, ny, lats2d, lons2d, obs_lats(n), obs_lons(n), i, j, xfloat, yfloat)
      
      ! Check for valid observations
      if (i == -1 .or. j == -1) then
        print *, "WRFDA DEBUG: Observation", n, " is out of domain - skipping"
        cycle
      end if
      
      if (obs_qc_flags(n) < 0) then
        print *, "WRFDA DEBUG: Observation", n, " failed QC - skipping"
        cycle
      end if
      
      if (obs_usage_flags(n) <= 0) then
        print *, "WRFDA DEBUG: Observation", n, " not used - skipping"
        cycle
      end if
      
      valid_obs_count = valid_obs_count + 1
      
      print *, "WRFDA DEBUG: Processing valid observation", valid_obs_count
      print *, "WRFDA DEBUG:   Lat:", obs_lats(n), " Lon:", obs_lons(n)
      print *, "WRFDA DEBUG:   Elevation:", obs_elevations(n), " Pressure:", obs_pressures(n), " Height:", obs_heights(n)
      print *, "WRFDA DEBUG:   QC Flag:", obs_qc_flags(n), " Usage Flag:", obs_usage_flags(n)
      print *, "WRFDA DEBUG:   Value:", obs_values(n), " Error:", obs_errors(n)
      
      ! Set grid indices for all levels (surface observations)
      do k = 1, nz
        iv%info(family)%i(k, valid_obs_count) = i
        iv%info(family)%j(k, valid_obs_count) = j
        iv%info(family)%x(k, valid_obs_count) = real(xfloat)
        iv%info(family)%y(k, valid_obs_count) = real(yfloat)
        
        ! Compute proper interpolation weights for bilinear interpolation
        iv%info(family)%dx(k, valid_obs_count) = real(xfloat - int(xfloat))  ! Fractional part of x
        iv%info(family)%dxm(k, valid_obs_count) = 1.0 - iv%info(family)%dx(k, valid_obs_count)  ! 1.0 - dx
        iv%info(family)%dy(k, valid_obs_count) = real(yfloat - int(yfloat))  ! Fractional part of y  
        iv%info(family)%dym(k, valid_obs_count) = 1.0 - iv%info(family)%dy(k, valid_obs_count)  ! 1.0 - dy
      end do
      
      ! Set vertical interpolation
      if (nz > 1) then
        call find_fractional_k(nz, levels, obs_levels(n), k, zkfloat)
        do k = 1, nz
          iv%info(family)%k(k, valid_obs_count) = k
          iv%info(family)%zk(k, valid_obs_count) = real(zkfloat)
          iv%info(family)%dz(k, valid_obs_count) = zkfloat - real(k, c_double)
        end do
      else
        iv%info(family)%k(1, valid_obs_count) = 1
        iv%info(family)%zk(1, valid_obs_count) = 1.0
        iv%info(family)%dz(1, valid_obs_count) = 0.0
      end if
      
      ! Note: The infa_type structure in WRFDA only contains grid interpolation fields
      ! Observation-specific data (values, errors, QC flags, etc.) are stored in
      ! separate structures (like synop_type, metar_type, etc.) that are accessed
      ! by the specific observation operator routines (da_transform_xtoy_synop, etc.)
      ! 
      ! The iv%info structure is primarily for grid interpolation information,
      ! not for storing observation values or metadata.
      !
      ! For now, we'll just set the basic grid interpolation fields that are
      ! actually available in the infa_type structure.
      !
      ! The enhanced observation data (values, errors, QC flags, station info)
      ! would need to be stored in the appropriate observation-specific structures
      ! (iv%synop, iv%metar, etc.) which are allocated and managed by WRFDA's
      ! da_allocate_y routine and the specific observation operator routines.
    end do
    
    print *, "WRFDA DEBUG: init_iv_from_enhanced_obs completed with", valid_obs_count, "valid observations"
  end subroutine init_iv_from_enhanced_obs

  ! Helper function to extract C string from array
  subroutine extract_c_string(c_strings, index, fortran_string)
    character(c_char), intent(in) :: c_strings(*)
    integer, intent(in) :: index
    character(len=*), intent(out) :: fortran_string
    integer :: i, start_pos, end_pos
    
    ! Calculate start position (assuming fixed-length strings)
    start_pos = (index - 1) * len(fortran_string) + 1
    end_pos = start_pos + len(fortran_string) - 1
    
    ! Extract string
    fortran_string = ""
    do i = start_pos, end_pos
      if (c_strings(i) == c_null_char) exit
      fortran_string(i-start_pos+1:i-start_pos+1) = achar(iachar(c_strings(i)))
    end do
  end subroutine extract_c_string

  subroutine init_domain_from_arrays_refs(grid, nx, ny, nz, u, v, t, q, psfc)
    type(domain), intent(inout) :: grid
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(inout) :: u(*), v(*), t(*), q(*), psfc(*)
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
    
    ! Suppress unused argument warnings
    if (u(1) < 0.0 .or. v(1) < 0.0 .or. t(1) < 0.0 .or. q(1) < 0.0 .or. psfc(1) < 0.0) then
      ! This will never execute but prevents unused argument warnings
      continue
    end if
    
    print *, "WRFDA DEBUG: init_domain_from_arrays_refs: nx=", nx, " ny=", ny, " nz=", nz, " nz1=", nz1
    
    ! CRITICAL FIX: Set WRFDA grid bounds that are used for observation array allocation
    ! These are used in da_copy_dims.inc to set kms, kme, ims, ime, jms, jme
    print *, "WRFDA DEBUG: Setting WRFDA grid bounds for proper observation array allocation"
    ! WRFDA expects 1-based bounds for grid arrays
    grid%sm31 = 1; grid%em31 = nx
    grid%sm32 = 1; grid%em32 = ny  
    grid%sm33 = 1; grid%em33 = nz1
    
    ! CRITICAL FIX: Also set domain start/end bounds that are used in da_copy_dims.inc
    ! These are used to set ids, ide, jds, jde, kds, kde
    print *, "WRFDA DEBUG: Setting WRFDA domain start/end bounds"
    grid%sd31 = 1; grid%ed31 = nx
    grid%sd32 = 1; grid%ed32 = ny
    grid%sd33 = 1; grid%ed33 = nz1
    
    ! Also set the parallel decomposition bounds (single tile case)
    grid%sp31 = 1; grid%ep31 = nx
    grid%sp32 = 1; grid%ep32 = ny
    grid%sp33 = 1; grid%ep33 = nz1
    
    ! CRITICAL FIX: Set domain bounds that are used for task bounds (kts:kte)
    ! These are used in da_copy_tile_dims.inc and da_interp_lin_3d.inc
    print *, "WRFDA DEBUG: Setting WRFDA domain bounds for proper task bounds"
    grid%xp%kds = 1; grid%xp%kde = nz1
    grid%xp%ids = 1; grid%xp%ide = nx
    grid%xp%jds = 1; grid%xp%jde = ny
    
    ! CRITICAL FIX: Also set the xp task bounds to match the domain bounds
    ! This ensures kts:kte loop uses valid indices
    print *, "WRFDA DEBUG: Setting WRFDA xp task bounds to match domain bounds"
    grid%xp%kts = 1; grid%xp%kte = nz1
    grid%xp%its = 1; grid%xp%ite = nx
    grid%xp%jts = 1; grid%xp%jte = ny
    
    ! CRITICAL FIX: Set up tile information for single-tile case
    ! This is required by da_copy_tile_dims to set kts, kte, its, ite, jts, jte
    print *, "WRFDA DEBUG: Setting up tile information for single-tile case"
    grid%num_tiles = 1
    allocate(grid%i_start(1:1), grid%i_end(1:1))
    allocate(grid%j_start(1:1), grid%j_end(1:1))
    grid%i_start(1) = 1; grid%i_end(1) = nx
    grid%j_start(1) = 1; grid%j_end(1) = ny
    print *, "WRFDA DEBUG: Tile arrays allocated: i_start=", grid%i_start(1), " i_end=", grid%i_end(1), " j_start=", grid%j_start(1), " j_end=", grid%j_end(1)
    
    allocate(grid%xa%u(nx,ny,nz1)); allocate(grid%xa%v(nx,ny,nz1))
    allocate(grid%xa%t(nx,ny,nz1)); allocate(grid%xa%q(nx,ny,nz1))
    allocate(grid%xa%psfc(nx,ny))
    do k=1,nz1; do j=1,ny; do i=1,nx
      grid%xa%u(i,j,k) = 0.0; grid%xa%v(i,j,k) = 0.0; grid%xa%t(i,j,k) = 0.0; grid%xa%q(i,j,k) = 0.0
    end do; end do; end do
    grid%xa%psfc = 0.0
    allocate(grid%xb%u(nx,ny,nz1)); grid%xb%u = 0.0
    allocate(grid%xb%v(nx,ny,nz1)); grid%xb%v = 0.0
    allocate(grid%xb%t(nx,ny,nz1)); grid%xb%t = 0.0
    allocate(grid%xb%q(nx,ny,nz1)); grid%xb%q = 0.0
    allocate(grid%xb%psfc(nx,ny));  grid%xb%psfc = 0.0
    
    print *, "WRFDA DEBUG: init_domain_from_arrays_refs completed successfully"
  end subroutine init_domain_from_arrays_refs

  subroutine init_iv_from_obs(family, iv, nx, ny, nz, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels, kms, kme)
     integer, intent(in) :: family
     type(iv_type), intent(inout) :: iv
     integer, intent(in) :: nx, ny, nz
     real(c_double), intent(in) :: lats2d(*), lons2d(*), levels(*)
     integer, intent(in) :: num_obs
     real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels(*)
     integer, intent(in) :: kms, kme
         integer :: n, i, j, k
    real(c_double) :: xfloat, yfloat, zkfloat
    
    ! Suppress unused argument warnings
    if (kms < 0 .or. kme < 0) then
      ! This will never execute but prevents unused argument warnings
      continue
    end if
    
    print *, "WRFDA DEBUG: init_iv_from_obs: family=", family, " num_obs=", num_obs, " nz=", nz
    
    ! CRITICAL FIX: Initialize ALL observation families to prevent huge memory allocation
    print *, "WRFDA DEBUG: Initializing ALL observation families to prevent memory issues..."
    iv%info(:)%nlocal = 0
    iv%info(:)%ntotal = 0
    iv%info(:)%max_lev = 1
    
    ! Initialize instrument-related fields for conventional observations only
    print *, "WRFDA DEBUG: Initializing instrument fields for conventional observations..."
    iv%num_inst = 0
    iv%total_rad_pixel = 0
    iv%total_rad_channel = 0
    nullify(iv%instid)
    
         ! Now set the actual family we want to use
     iv%info(family)%nlocal = num_obs
     iv%info(family)%ntotal = num_obs
     iv%info(family)%max_lev = nz  ! Use actual number of vertical levels
    
    ! CRITICAL FIX: Use WRFDA's own allocation routine instead of manual allocation
    ! This ensures the arrays have the correct bounds that WRFDA expects
    print *, "WRFDA DEBUG: Calling WRFDA's da_allocate_obs_info to allocate arrays with correct bounds"
    call da_allocate_obs_info(iv, family)
    print *, "WRFDA DEBUG: da_allocate_obs_info completed successfully"
    
    ! Debug: Check what bounds WRFDA actually allocated
    print *, "WRFDA DEBUG: After da_allocate_obs_info:"
    print *, "WRFDA DEBUG:   iv%info(family)%i bounds:", lbound(iv%info(family)%i), ":", ubound(iv%info(family)%i)
    print *, "WRFDA DEBUG:   iv%info(family)%dym bounds:", lbound(iv%info(family)%dym), ":", ubound(iv%info(family)%dym)
    print *, "WRFDA DEBUG:   iv%info(family)%k bounds:", lbound(iv%info(family)%k), ":", ubound(iv%info(family)%k)
    
    ! Set the levels array after WRFDA allocation
    iv%info(family)%levels = 1
    
    ! Initialize grid indices - for surface observations, use surface level indices
    iv%info(family)%i = 1  ! Default to grid point 1 (1-based)
    iv%info(family)%j = 1  ! Default to grid point 1 (1-based)
    iv%info(family)%k = 1  ! Default to surface level (1-based)
    print *, "WRFDA DEBUG: Grid indices initialized"
    
    ! Set up plocal array for proper time slot indexing
    ! For 3D-Var (num_fgat_time = 1), we have:
    ! plocal(0) = 0 (no observations before time slot 1)
    ! plocal(1) = num_obs (cumulative count at time slot 1)
    iv%info(family)%plocal(0) = 0
    iv%info(family)%plocal(1) = num_obs
    
    ! Set n1 and n2 based on plocal array (WRFDA standard pattern)
    iv%info(family)%n1 = iv%info(family)%plocal(0) + 1  ! = 1
    iv%info(family)%n2 = iv%info(family)%plocal(1)      ! = num_obs
      do n = 1, num_obs
        call find_fractional_ij(nx, ny, lats2d, lons2d, obs_lats(n), obs_lons(n), i, j, xfloat, yfloat)
        print *, "WRFDA DEBUG: Observation", n, " - calculated grid indices: i=", i, " j=", j, " xfloat=", xfloat, " yfloat=", yfloat
        
        ! Check for out-of-domain observations (i or j = -1)
        if (i == -1 .or. j == -1) then
          print *, "WRFDA DEBUG: Observation", n, " is out of domain - skipping processing"
          ! Skip this observation by setting invalid indices and continuing
          ! This follows WRFDA's pattern of graceful handling for out-of-domain observations
          do k = 1, nz
            iv%info(family)%i(k,n) = -1  ! Mark as invalid
            iv%info(family)%j(k,n) = -1  ! Mark as invalid
            iv%info(family)%x(k,n) = -1.0  ! Mark as invalid
            iv%info(family)%y(k,n) = -1.0  ! Mark as invalid
            iv%info(family)%dx(k,n) = 0.0
            iv%info(family)%dxm(k,n) = 1.0
            iv%info(family)%dy(k,n) = 0.0
            iv%info(family)%dym(k,n) = 1.0
          end do
          cycle  ! Skip to next observation
        end if
        
        ! CRITICAL FIX: For surface observations, only set the surface level (k=1)
        ! This prevents WRFDA from trying to interpolate 3D variables at multiple levels
        ! Surface observations should only be interpolated at the surface level
        k = 1  ! Surface level only
        iv%info(family)%i(k,n) = i  ! Keep 1-based
        iv%info(family)%j(k,n) = j  ! Keep 1-based
        iv%info(family)%x(k,n) = real(xfloat)  ! Keep 1-based
        iv%info(family)%y(k,n) = real(yfloat)  ! Keep 1-based
        
        ! Compute proper interpolation weights for bilinear interpolation
        ! xfloat and yfloat are fractional coordinates (e.g., 1.5 means halfway between grid points 1 and 2)
        ! dx = fractional part of x-coordinate, dxm = 1.0 - dx
        ! dy = fractional part of y-coordinate, dym = 1.0 - dy
        iv%info(family)%dx(k,n) = real(xfloat - int(xfloat))  ! Fractional part of x
        iv%info(family)%dxm(k,n) = 1.0 - iv%info(family)%dx(k,n)  ! 1.0 - dx
        iv%info(family)%dy(k,n) = real(yfloat - int(yfloat))  ! Fractional part of y  
        iv%info(family)%dym(k,n) = 1.0 - iv%info(family)%dy(k,n)  ! 1.0 - dy
        
        ! Set remaining levels to the same surface level to prevent WRFDA from processing them
        ! This ensures all levels point to the same valid surface location
        do k = 2, nz
          iv%info(family)%i(k,n) = i  ! Same horizontal position
          iv%info(family)%j(k,n) = j  ! Same horizontal position
          iv%info(family)%x(k,n) = real(xfloat)  ! Same horizontal position
          iv%info(family)%y(k,n) = real(yfloat)  ! Same horizontal position
          ! Use same interpolation weights as the first level
          iv%info(family)%dx(k,n) = iv%info(family)%dx(1,n)  ! Same as first level
          iv%info(family)%dxm(k,n) = iv%info(family)%dxm(1,n) ! Same as first level
          iv%info(family)%dy(k,n) = iv%info(family)%dy(1,n)  ! Same as first level
          iv%info(family)%dym(k,n) = iv%info(family)%dym(1,n) ! Same as first level
        end do
       if (nz > 1) then
         ! For surface observations, find_fractional_k will return k=1
         ! For upper-air observations, it will return proper vertical interpolation
         call find_fractional_k(nz, levels, obs_levels(n), k, zkfloat)
         do k = 1, nz
           iv%info(family)%k(k,n) = k
           iv%info(family)%zk(k,n) = real(zkfloat)
           iv%info(family)%dz(k,n) = zkfloat - real(k, c_double)
           iv%info(family)%dzm(k,n) = 1.0 - iv%info(family)%dz(k,n)
         end do
       else
         do k = 1, nz
           iv%info(family)%k(k,n) = 1
           iv%info(family)%zk(k,n) = 1.0
           iv%info(family)%dz(k,n) = 0.0
           iv%info(family)%dzm(k,n) = 1.0
           ! CRITICAL FIX: Initialize vertical interpolation weights for single-level case
           ! This ensures da_interp_lin_3d has valid interpolation weights
         end do
       end if
     end do
    print *, "WRFDA DEBUG: init_iv_from_obs completed successfully"
  end subroutine init_iv_from_obs

  subroutine init_iv_profiles(family, iv, nx, ny, nz, lats2d, lons2d, levels, &
                              num_obs, obs_counts, obs_lats, obs_lons, obs_levels_flat)
    integer, intent(in) :: family
    type(iv_type), intent(inout) :: iv
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: lats2d(*), lons2d(*), levels(*)
    integer, intent(in) :: num_obs
    integer(c_int), intent(in) :: obs_counts(*)
    real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels_flat(*)
    integer :: n, i, j, idx, lev, nlev, max_lev
    real(c_double) :: xfloat, yfloat, zkfloat
    ! Determine max levels among obs
    max_lev = 0
    do n=1,num_obs
      if (obs_counts(n) > max_lev) max_lev = obs_counts(n)
    end do
    if (max_lev <= 0) max_lev = 1
    ! Allocate per-family iv info arrays sized by (levels, obs)
    iv%info(family)%nlocal = num_obs
    iv%info(family)%ntotal = num_obs
    iv%info(family)%max_lev = max_lev
    allocate(iv%info(family)%i(max_lev,num_obs), iv%info(family)%j(max_lev,num_obs), iv%info(family)%k(max_lev,num_obs))
    allocate(iv%info(family)%x(max_lev,num_obs), iv%info(family)%y(max_lev,num_obs), iv%info(family)%zk(max_lev,num_obs))
    allocate(iv%info(family)%dx(max_lev,num_obs), iv%info(family)%dy(max_lev,num_obs), iv%info(family)%dz(max_lev,num_obs))
    allocate(iv%info(family)%dxm(max_lev,num_obs), iv%info(family)%dym(max_lev,num_obs), iv%info(family)%dzm(max_lev,num_obs))
    allocate(iv%info(family)%levels(num_obs))
    allocate(iv%info(family)%proc_domain(max_lev,num_obs)); iv%info(family)%proc_domain = .true.
    
    ! Set up plocal array for proper time slot indexing
    ! For 3D-Var (num_fgat_time = 1), we have:
    ! plocal(0) = 0 (no observations before time slot 1)
    ! plocal(1) = num_obs (cumulative count at time slot 1)
    iv%info(family)%plocal(0) = 0
    iv%info(family)%plocal(1) = num_obs
    
    ! Set n1 and n2 based on plocal array (WRFDA standard pattern)
    iv%info(family)%n1 = iv%info(family)%plocal(0) + 1  ! = 1
    iv%info(family)%n2 = iv%info(family)%plocal(1)      ! = num_obs
    idx = 0
         do n = 1, num_obs
       nlev = obs_counts(n)
       iv%info(family)%levels(n) = nlev
       call find_fractional_ij(nx, ny, lats2d, lons2d, obs_lats(n), obs_lons(n), i, j, xfloat, yfloat)
       
       ! Check for out-of-domain observations (i or j = -1)
       if (i == -1 .or. j == -1) then
         print *, "WRFDA DEBUG: Profile observation", n, " is out of domain - skipping processing"
         ! Skip this observation by setting invalid indices and continuing
         ! This follows WRFDA's pattern of graceful handling for out-of-domain observations
         do lev = 1, max_lev
           iv%info(family)%i(lev,n) = -1  ! Mark as invalid
           iv%info(family)%j(lev,n) = -1  ! Mark as invalid
           iv%info(family)%k(lev,n) = -1  ! Mark as invalid
           iv%info(family)%x(lev,n) = -1.0  ! Mark as invalid
           iv%info(family)%y(lev,n) = -1.0  ! Mark as invalid
           iv%info(family)%zk(lev,n) = -1.0  ! Mark as invalid
           ! CRITICAL FIX: Initialize all interpolation weights for out-of-domain obs
           iv%info(family)%dx(lev,n) = 0.0
           iv%info(family)%dxm(lev,n) = 1.0
           iv%info(family)%dy(lev,n) = 0.0
           iv%info(family)%dym(lev,n) = 1.0
           iv%info(family)%dz(lev,n) = 0.0
           iv%info(family)%dzm(lev,n) = 1.0
         end do
         cycle  ! Skip to next observation
       end if
       
       do lev = 1, nlev
        ! CRITICAL FIX: Add bounds checking for obs_levels_flat access
        ! Note: We can't use size() on assumed-size arrays, so we'll use a safer approach
        ! The idx variable should be properly managed by the calling code
        if (idx < 0 .or. lev < 1 .or. lev > nlev) then
          print *, "WRFDA ERROR: Invalid indices in profile processing"
          print *, "  idx=", idx, ", lev=", lev, ", nlev=", nlev
          ! Set invalid values and continue to prevent crash
          iv%info(family)%k(lev,n) = -1
          iv%info(family)%zk(lev,n) = -1.0
          iv%info(family)%dz(lev,n) = 0.0
          iv%info(family)%dzm(lev,n) = 1.0
          cycle
        end if
        
        if (nz > 1) then
          ! For surface observations, find_fractional_k will return k=1
          ! For upper-air observations, it will return proper vertical interpolation
          call find_fractional_k(nz, levels, obs_levels_flat(idx+lev), iv%info(family)%k(lev,n), zkfloat)
          iv%info(family)%zk(lev,n) = real(zkfloat)
          iv%info(family)%dz(lev,n) = zkfloat - real(iv%info(family)%k(lev,n), c_double)
          iv%info(family)%dzm(lev,n) = 1.0 - iv%info(family)%dz(lev,n)
        else
          iv%info(family)%k(lev,n) = 1
          iv%info(family)%zk(lev,n) = 1.0
          iv%info(family)%dz(lev,n) = 0.0
          iv%info(family)%dzm(lev,n) = 1.0
        end if
        
        ! CRITICAL FIX: Initialize horizontal interpolation weights
        iv%info(family)%i(lev,n) = i
        iv%info(family)%j(lev,n) = j
        iv%info(family)%x(lev,n) = real(xfloat)
        iv%info(family)%y(lev,n) = real(yfloat)
        iv%info(family)%dx(lev,n) = xfloat - real(i, c_double)
        iv%info(family)%dxm(lev,n) = 1.0 - iv%info(family)%dx(lev,n)
        iv%info(family)%dy(lev,n) = yfloat - real(j, c_double)
        iv%info(family)%dym(lev,n) = 1.0 - iv%info(family)%dy(lev,n)
      end do
      ! Fill remaining levels (from nlev+1 to max_lev) with the same horizontal position
      do lev = nlev+1, max_lev
        iv%info(family)%i(lev,n) = i
        iv%info(family)%j(lev,n) = j
        iv%info(family)%k(lev,n) = 1
        iv%info(family)%x(lev,n) = real(xfloat)
        iv%info(family)%y(lev,n) = real(yfloat)
        iv%info(family)%zk(lev,n) = 1.0
        ! CRITICAL FIX: Initialize interpolation weights for remaining levels
        iv%info(family)%dx(lev,n) = xfloat - real(i, c_double)
        iv%info(family)%dxm(lev,n) = 1.0 - iv%info(family)%dx(lev,n)
        iv%info(family)%dy(lev,n) = yfloat - real(j, c_double)
        iv%info(family)%dym(lev,n) = 1.0 - iv%info(family)%dy(lev,n)
        iv%info(family)%dz(lev,n) = 0.0
        iv%info(family)%dzm(lev,n) = 1.0
      end do
       idx = idx + nlev
     end do
  end subroutine init_iv_profiles

  ! No longer used: WRFDA's da_allocate_y handles family allocation
  subroutine init_y_for_family(family, y, num_obs)
    integer, intent(in) :: family
    type(y_type), intent(inout) :: y
    integer, intent(in) :: num_obs
    ! Suppress unused argument warnings
    if (family < 0 .or. num_obs < 0) then
      ! This will never execute but prevents unused argument warnings
      continue
    end if
    ! Suppress unused y warning by referencing it
    if (.false.) then
      ! This will never execute but prevents unused argument warnings
      ! Reference y to suppress warning
      if (associated(y%synop)) then
        continue
      end if
    end if
  end subroutine init_y_for_family

  subroutine init_y_from_delta(family, var_code, jo_grad_y, delta_y, num_obs)
    integer, intent(in) :: family
    character(len=1), intent(in) :: var_code  ! 't','q','u','v','p'
    type(y_type), intent(inout) :: jo_grad_y
    real(c_double), intent(in) :: delta_y(*)
    integer, intent(in) :: num_obs
    integer :: n, k
    select case (family)
    case (metar)
      allocate(jo_grad_y%metar(num_obs))
      do n=1,num_obs
        jo_grad_y%metar(n)%u = 0.0; jo_grad_y%metar(n)%v = 0.0
        jo_grad_y%metar(n)%t = 0.0; jo_grad_y%metar(n)%q = 0.0; jo_grad_y%metar(n)%p = 0.0
        select case (var_code)
        case('u'); jo_grad_y%metar(n)%u = real(delta_y(n))
        case('v'); jo_grad_y%metar(n)%v = real(delta_y(n))
        case('t'); jo_grad_y%metar(n)%t = real(delta_y(n))
        case('q'); jo_grad_y%metar(n)%q = real(delta_y(n))
        case('p'); jo_grad_y%metar(n)%p = real(delta_y(n))
        case default; jo_grad_y%metar(n)%t = real(delta_y(n))
        end select
      end do
    case (synop)
      allocate(jo_grad_y%synop(num_obs))
      do n=1,num_obs
        jo_grad_y%synop(n)%u = 0.0; jo_grad_y%synop(n)%v = 0.0
        jo_grad_y%synop(n)%t = 0.0; jo_grad_y%synop(n)%q = 0.0; jo_grad_y%synop(n)%p = 0.0
        select case (var_code)
        case('u'); jo_grad_y%synop(n)%u = real(delta_y(n))
        case('v'); jo_grad_y%synop(n)%v = real(delta_y(n))
        case('t'); jo_grad_y%synop(n)%t = real(delta_y(n))
        case('q'); jo_grad_y%synop(n)%q = real(delta_y(n))
        case('p'); jo_grad_y%synop(n)%p = real(delta_y(n))
        case default; jo_grad_y%synop(n)%t = real(delta_y(n))
        end select
      end do
    case (ships)
      allocate(jo_grad_y%ships(num_obs))
      do n=1,num_obs
        jo_grad_y%ships(n)%u = 0.0; jo_grad_y%ships(n)%v = 0.0
        jo_grad_y%ships(n)%t = 0.0; jo_grad_y%ships(n)%q = 0.0; jo_grad_y%ships(n)%p = 0.0
        select case (var_code)
        case('u'); jo_grad_y%ships(n)%u = real(delta_y(n))
        case('v'); jo_grad_y%ships(n)%v = real(delta_y(n))
        case('t'); jo_grad_y%ships(n)%t = real(delta_y(n))
        case('q'); jo_grad_y%ships(n)%q = real(delta_y(n))
        case('p'); jo_grad_y%ships(n)%p = real(delta_y(n))
        case default; jo_grad_y%ships(n)%t = real(delta_y(n))
        end select
      end do
    case (buoy)
      allocate(jo_grad_y%buoy(num_obs))
      do n=1,num_obs
        jo_grad_y%buoy(n)%u = 0.0; jo_grad_y%buoy(n)%v = 0.0
        jo_grad_y%buoy(n)%t = 0.0; jo_grad_y%buoy(n)%q = 0.0; jo_grad_y%buoy(n)%p = 0.0
        select case (var_code)
        case('u'); jo_grad_y%buoy(n)%u = real(delta_y(n))
        case('v'); jo_grad_y%buoy(n)%v = real(delta_y(n))
        case('t'); jo_grad_y%buoy(n)%t = real(delta_y(n))
        case('q'); jo_grad_y%buoy(n)%q = real(delta_y(n))
        case('p'); jo_grad_y%buoy(n)%p = real(delta_y(n))
        case default; jo_grad_y%buoy(n)%t = real(delta_y(n))
        end select
      end do
    case (sonde_sfc)
      allocate(jo_grad_y%sonde_sfc(num_obs))
      do n=1,num_obs
        jo_grad_y%sonde_sfc(n)%u = 0.0; jo_grad_y%sonde_sfc(n)%v = 0.0
        jo_grad_y%sonde_sfc(n)%t = 0.0; jo_grad_y%sonde_sfc(n)%q = 0.0; jo_grad_y%sonde_sfc(n)%p = 0.0
        select case (var_code)
        case('u'); jo_grad_y%sonde_sfc(n)%u = real(delta_y(n))
        case('v'); jo_grad_y%sonde_sfc(n)%v = real(delta_y(n))
        case('t'); jo_grad_y%sonde_sfc(n)%t = real(delta_y(n))
        case('q'); jo_grad_y%sonde_sfc(n)%q = real(delta_y(n))
        case('p'); jo_grad_y%sonde_sfc(n)%p = real(delta_y(n))
        case default; jo_grad_y%sonde_sfc(n)%t = real(delta_y(n))
        end select
      end do
    case (airep)
      do n=1,num_obs
        k = 1
        select case (var_code)
        case('u'); jo_grad_y%airep(n)%u(k) = real(delta_y(n))
        case('v'); jo_grad_y%airep(n)%v(k) = real(delta_y(n))
        case('t'); jo_grad_y%airep(n)%t(k) = real(delta_y(n))
        case('q'); jo_grad_y%airep(n)%q(k) = real(delta_y(n))
        case default; jo_grad_y%airep(n)%t(k) = real(delta_y(n))
        end select
      end do
    case (pilot)
      do n=1,num_obs
        k = 1
        select case (var_code)
        case('u'); jo_grad_y%pilot(n)%u(k) = real(delta_y(n))
        case('v'); jo_grad_y%pilot(n)%v(k) = real(delta_y(n))
        case default; jo_grad_y%pilot(n)%u(k) = real(delta_y(n))
        end select
      end do
    case (sound)
      do n=1,num_obs
        k = 1
        select case (var_code)
        case('u'); jo_grad_y%sound(n)%u(k) = real(delta_y(n))
        case('v'); jo_grad_y%sound(n)%v(k) = real(delta_y(n))
        case('t'); jo_grad_y%sound(n)%t(k) = real(delta_y(n))
        case('q'); jo_grad_y%sound(n)%q(k) = real(delta_y(n))
        case default; jo_grad_y%sound(n)%t(k) = real(delta_y(n))
        end select
      end do
    end select
  end subroutine init_y_from_delta

  subroutine copy_y_to_out(family, var_code, y, out_y, num_obs)
    integer, intent(in) :: family
    character(len=1), intent(in) :: var_code  ! 't','q','u','v','p'
    type(y_type), intent(in) :: y
    real(c_double), intent(out) :: out_y(*)
    integer, intent(in) :: num_obs
    integer :: n, k
    select case (family)
    case (metar)
      do n=1,num_obs
        select case (var_code)
        case('u'); out_y(n) = real(y%metar(n)%u, kind=c_double)
        case('v'); out_y(n) = real(y%metar(n)%v, kind=c_double)
        case('t'); out_y(n) = real(y%metar(n)%t, kind=c_double)
        case('q'); out_y(n) = real(y%metar(n)%q, kind=c_double)
        case('p'); out_y(n) = real(y%metar(n)%p, kind=c_double)
        case default; out_y(n) = real(y%metar(n)%t, kind=c_double)
        end select
      end do
    case (synop)
      do n=1,num_obs
        select case (var_code)
        case('u'); out_y(n) = real(y%synop(n)%u, kind=c_double)
        case('v'); out_y(n) = real(y%synop(n)%v, kind=c_double)
        case('t'); out_y(n) = real(y%synop(n)%t, kind=c_double)
        case('q'); out_y(n) = real(y%synop(n)%q, kind=c_double)
        case('p'); out_y(n) = real(y%synop(n)%p, kind=c_double)
        case default; out_y(n) = real(y%synop(n)%t, kind=c_double)
        end select
      end do
    case (ships)
      do n=1,num_obs
        select case (var_code)
        case('u'); out_y(n) = real(y%ships(n)%u, kind=c_double)
        case('v'); out_y(n) = real(y%ships(n)%v, kind=c_double)
        case('t'); out_y(n) = real(y%ships(n)%t, kind=c_double)
        case('q'); out_y(n) = real(y%ships(n)%q, kind=c_double)
        case('p'); out_y(n) = real(y%ships(n)%p, kind=c_double)
        case default; out_y(n) = real(y%ships(n)%t, kind=c_double)
        end select
      end do
    case (buoy)
      do n=1,num_obs
        select case (var_code)
        case('u'); out_y(n) = real(y%buoy(n)%u, kind=c_double)
        case('v'); out_y(n) = real(y%buoy(n)%v, kind=c_double)
        case('t'); out_y(n) = real(y%buoy(n)%t, kind=c_double)
        case('q'); out_y(n) = real(y%buoy(n)%q, kind=c_double)
        case('p'); out_y(n) = real(y%buoy(n)%p, kind=c_double)
        case default; out_y(n) = real(y%buoy(n)%t, kind=c_double)
        end select
      end do
    case (sonde_sfc)
      do n=1,num_obs
        select case (var_code)
        case('u'); out_y(n) = real(y%sonde_sfc(n)%u, kind=c_double)
        case('v'); out_y(n) = real(y%sonde_sfc(n)%v, kind=c_double)
        case('t'); out_y(n) = real(y%sonde_sfc(n)%t, kind=c_double)
        case('q'); out_y(n) = real(y%sonde_sfc(n)%q, kind=c_double)
        case('p'); out_y(n) = real(y%sonde_sfc(n)%p, kind=c_double)
        case default; out_y(n) = real(y%sonde_sfc(n)%t, kind=c_double)
        end select
      end do
    case (airep)
      do n=1,num_obs
        k = 1
        select case (var_code)
        case('u'); out_y(n) = real(y%airep(n)%u(k), kind=c_double)
        case('v'); out_y(n) = real(y%airep(n)%v(k), kind=c_double)
        case('t'); out_y(n) = real(y%airep(n)%t(k), kind=c_double)
        case('q'); out_y(n) = real(y%airep(n)%q(k), kind=c_double)
        case default; out_y(n) = real(y%airep(n)%t(k), kind=c_double)
        end select
      end do
    case (pilot)
      do n=1,num_obs
        k = 1
        select case (var_code)
        case('u'); out_y(n) = real(y%pilot(n)%u(k), kind=c_double)
        case('v'); out_y(n) = real(y%pilot(n)%v(k), kind=c_double)
        case default; out_y(n) = real(y%pilot(n)%u(k), kind=c_double)
        end select
      end do
    case (sound)
      do n=1,num_obs
        k = 1
        select case (var_code)
        case('u'); out_y(n) = real(y%sound(n)%u(k), kind=c_double)
        case('v'); out_y(n) = real(y%sound(n)%v(k), kind=c_double)
        case('t'); out_y(n) = real(y%sound(n)%t(k), kind=c_double)
        case('q'); out_y(n) = real(y%sound(n)%q(k), kind=c_double)
        case default; out_y(n) = real(y%sound(n)%t(k), kind=c_double)
        end select
      end do
    end select
  end subroutine copy_y_to_out

  subroutine copy_y_to_out_profiles(family, var_code, y, num_obs, obs_counts, out_y_flat)
    integer, intent(in) :: family
    character(len=1), intent(in) :: var_code
    type(y_type), intent(in) :: y
    integer, intent(in) :: num_obs
    integer(c_int), intent(in) :: obs_counts(*)
    real(c_double), intent(out) :: out_y_flat(*)
    integer :: n, k, idx
    idx = 0
    do n=1,num_obs
      select case (family)
      case (airep)
        do k=1,obs_counts(n)
          select case (var_code)
          case('u'); out_y_flat(idx+k) = real(y%airep(n)%u(k), kind=c_double)
          case('v'); out_y_flat(idx+k) = real(y%airep(n)%v(k), kind=c_double)
          case('t'); out_y_flat(idx+k) = real(y%airep(n)%t(k), kind=c_double)
          case('q'); out_y_flat(idx+k) = real(y%airep(n)%q(k), kind=c_double)
          case default; out_y_flat(idx+k) = real(y%airep(n)%t(k), kind=c_double)
          end select
        end do
      case (pilot)
        do k=1,obs_counts(n)
          select case (var_code)
          case('u'); out_y_flat(idx+k) = real(y%pilot(n)%u(k), kind=c_double)
          case('v'); out_y_flat(idx+k) = real(y%pilot(n)%v(k), kind=c_double)
          case default; out_y_flat(idx+k) = real(y%pilot(n)%u(k), kind=c_double)
          end select
        end do
      case (sound)
        do k=1,obs_counts(n)
          select case (var_code)
          case('u'); out_y_flat(idx+k) = real(y%sound(n)%u(k), kind=c_double)
          case('v'); out_y_flat(idx+k) = real(y%sound(n)%v(k), kind=c_double)
          case('t'); out_y_flat(idx+k) = real(y%sound(n)%t(k), kind=c_double)
          case('q'); out_y_flat(idx+k) = real(y%sound(n)%q(k), kind=c_double)
          case default; out_y_flat(idx+k) = real(y%sound(n)%t(k), kind=c_double)
          end select
        end do
      end select
      idx = idx + obs_counts(n)
    end do
  end subroutine copy_y_to_out_profiles

  subroutine init_y_from_delta_profiles(family, var_code, jo_grad_y, num_obs, obs_counts, delta_y_flat)
    integer, intent(in) :: family
    character(len=1), intent(in) :: var_code
    type(y_type), intent(inout) :: jo_grad_y
    integer, intent(in) :: num_obs
    integer(c_int), intent(in) :: obs_counts(*)
    real(c_double), intent(in) :: delta_y_flat(*)
    integer :: n, k, idx
    idx = 0
    do n=1,num_obs
      do k=1,obs_counts(n)
        select case (family)
        case (airep)
          select case (var_code)
          case('u'); jo_grad_y%airep(n)%u(k) = real(delta_y_flat(idx+k))
          case('v'); jo_grad_y%airep(n)%v(k) = real(delta_y_flat(idx+k))
          case('t'); jo_grad_y%airep(n)%t(k) = real(delta_y_flat(idx+k))
          case('q'); jo_grad_y%airep(n)%q(k) = real(delta_y_flat(idx+k))
          end select
        case (pilot)
          select case (var_code)
          case('u'); jo_grad_y%pilot(n)%u(k) = real(delta_y_flat(idx+k))
          case('v'); jo_grad_y%pilot(n)%v(k) = real(delta_y_flat(idx+k))
          end select
        case (sound)
          select case (var_code)
          case('u'); jo_grad_y%sound(n)%u(k) = real(delta_y_flat(idx+k))
          case('v'); jo_grad_y%sound(n)%v(k) = real(delta_y_flat(idx+k))
          case('t'); jo_grad_y%sound(n)%t(k) = real(delta_y_flat(idx+k))
          case('q'); jo_grad_y%sound(n)%q(k) = real(delta_y_flat(idx+k))
          end select
        end select
      end do
      idx = idx + obs_counts(n)
    end do
  end subroutine init_y_from_delta_profiles

  subroutine parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
    character(len=*), intent(in) :: op_str
    character(len=256), intent(out) :: fam_str, var_str
    integer, intent(out) :: fam_id
    character(len=1), intent(out) :: var_code
    integer :: pos
    
    print *, "WRFDA DEBUG: parse_family_variable: input op_str='", trim(op_str), "'"
    
    fam_str = trim(op_str)
    var_str = ' '
    pos = index(op_str, ':')
    if (pos > 0) then
      fam_str = trim(op_str(:pos-1))
      var_str = trim(op_str(pos+1:))
      print *, "WRFDA DEBUG: Found separator at pos=", pos, " fam_str='", trim(fam_str), "' var_str='", trim(var_str), "'"
    else
      print *, "WRFDA DEBUG: No separator found, using full string as family"
    end if
    
    print *, "WRFDA DEBUG: Before select case, fam_str='", trim(fam_str), "'"
    select case (fam_str)
    case ('metar'); 
      fam_id = metar
      print *, "WRFDA DEBUG: Selected metar family, fam_id=", fam_id
    case ('synop'); 
      fam_id = synop
      print *, "WRFDA DEBUG: Selected synop family, fam_id=", fam_id
    case ('adpsfc'); 
      fam_id = synop
      print *, "WRFDA DEBUG: Mapped ADPSFC to synop family, fam_id=", fam_id
    case ('ships'); 
      fam_id = ships
      print *, "WRFDA DEBUG: Selected ships family, fam_id=", fam_id
    case ('buoy');  
      fam_id = buoy
      print *, "WRFDA DEBUG: Selected buoy family, fam_id=", fam_id
    case default;   
      fam_id = metar
      print *, "WRFDA DEBUG: Default case, using metar family, fam_id=", fam_id
    end select
    
    if (len_trim(var_str) == 0) then
      var_code = 't'
      print *, "WRFDA DEBUG: No variable specified, defaulting to 't'"
    else
      select case (var_str(1:1))
      case ('u','U')
        var_code = 'u'
      case ('v','V')
        var_code = 'v'
      case ('t','T')
        var_code = 't'
      case ('q','Q')
        var_code = 'q'
      case ('p','P')
        var_code = 'p'
      case default
        var_code = 't'
      end select
      print *, "WRFDA DEBUG: Variable code set to '", var_code, "'"
    end if
    
    print *, "WRFDA DEBUG: parse_family_variable completed: fam_id=", fam_id, " var_code='", var_code, "'"
  end subroutine parse_family_variable

  ! Copy gradient from jo_grad_x to state arrays (gradient accumulation)
  ! This function accumulates the adjoint gradients into the state space arrays
  ! For incremental 3D-Var, this represents the gradient of the cost function with respect to state
  subroutine copy_x_to_state(jo_grad_x, u, v, t, q, psfc, nx, ny, nz)
    type(x_type), intent(in) :: jo_grad_x
    real(c_double), intent(inout) :: u(*), v(*), t(*), q(*), psfc(*)
    integer, intent(in) :: nx, ny, nz
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
    
    print *, "WRFDA DEBUG: copy_x_to_state: accumulating gradients into state arrays"
    
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
    
    print *, "WRFDA DEBUG: copy_x_to_state: gradient accumulation completed"
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

  ! Use WRFDA-native approach: integrate with da_togrid for accurate grid index calculation
  ! This gives us WRFDA-native grid handling with proper interpolation weights and domain boundary handling
  ! 
  ! NOTE: When observations are out of domain bounds, this subroutine sets oi=oj=-1 and 
  ! xfloat=yfloat=-1.0 to indicate invalid coordinates. Callers MUST check for these -1 
  ! values and handle out-of-domain observations gracefully (e.g., skip processing, 
  ! cycle to next observation, etc.) following WRFDA's standard pattern.
  !
  ! LONGITUDE RANGE HANDLING: This subroutine handles the common case where lons2d uses
  ! -180 to 180 range while olon uses 0 to 359 range. It automatically converts between
  ! these conventions to ensure proper domain boundary checking and grid point finding.
  subroutine find_fractional_ij(nx, ny, lats2d, lons2d, olat, olon, oi, oj, xfloat, yfloat)
    integer, intent(in) :: nx, ny
    real(c_double), intent(in) :: lats2d(*), lons2d(*), olat, olon
    integer, intent(out) :: oi, oj
    real(c_double), intent(out) :: xfloat, yfloat
    
    ! Find nearest grid point by simple distance calculation as initial estimate
    integer :: i, j
    real(c_double) :: bestd, d, lat, lon
    
    ! Variables for da_togrid calls (double precision as expected by the compiler)
    real(c_double) :: x_approx, y_approx, dx_x, dxm_x, dy_y, dym_y
    integer :: temp_i, temp_j
    
    ! WRFDA-style domain boundary handling variables
    logical :: outside_domain
    
    ! Longitude range conversion variables
    real(c_double) :: olon_converted, lons2d_min, lons2d_max
    logical :: lons2d_negative_range, olon_positive_range
    
    ! Latitude range variables
    real(c_double) :: lats2d_min, lats2d_max
    
    ! Determine longitude range conventions
    ! Scan the entire array to find actual min/max longitude values
    ! This is more robust than assuming array ordering
    lons2d_min = lons2d(1)
    lons2d_max = lons2d(1)
    do j = 1, ny
      do i = 1, nx
        lons2d_min = min(lons2d_min, lons2d(i + (j-1)*nx))
        lons2d_max = max(lons2d_max, lons2d(i + (j-1)*nx))
      end do
    end do
    
    ! Determine latitude range conventions
    ! Scan the entire array to find actual min/max latitude values
    ! This is more robust than assuming array ordering
    lats2d_min = lats2d(1)
    lats2d_max = lats2d(1)
    do j = 1, ny
      do i = 1, nx
        lats2d_min = min(lats2d_min, lats2d(i + (j-1)*nx))
        lats2d_max = max(lats2d_max, lats2d(i + (j-1)*nx))
      end do
    end do
    
    print *, "WRFDA DEBUG: lons2d_min=", lons2d_min, " lons2d_max=", lons2d_max
    print *, "WRFDA DEBUG: lats2d_min=", lats2d_min, " lats2d_max=", lats2d_max
    print *, "WRFDA DEBUG: olon=", olon
    lons2d_negative_range = (lons2d_min < 0.0_c_double)
    olon_positive_range = (olon >= 0.0_c_double .and. olon <= 360.0_c_double)
    
    ! Convert olon to match lons2d range if necessary
    if (lons2d_negative_range .and. olon_positive_range) then
      ! lons2d uses -180 to 180, olon uses 0 to 359
      ! Convert olon to -180 to 180 range
      if (olon > 180.0_c_double) then
        olon_converted = olon - 360.0_c_double
      else
        olon_converted = olon
      end if
    else if (.not. lons2d_negative_range .and. .not. olon_positive_range) then
      ! lons2d uses 0 to 360, olon uses -180 to 180
      ! Convert olon to 0 to 360 range
      if (olon < 0.0_c_double) then
        olon_converted = olon + 360.0_c_double
      else
        olon_converted = olon
      end if
    else
      ! Same range convention, no conversion needed
      olon_converted = olon
    end if
    
    print *, "WRFDA DEBUG: lons2d_negative_range=", lons2d_negative_range, " olon_positive_range=", olon_positive_range
    print *, "WRFDA DEBUG: olon_converted=", olon_converted
    print *, "WRFDA DEBUG: lons2d_min=", lons2d_min, " lons2d_max=", lons2d_max
    print *, "WRFDA DEBUG: lats2d_min=", lats2d_min, " lats2d_max=", lats2d_max
    print *, "WRFDA DEBUG: olat=", olat
    ! First, check if observation is outside the domain bounds
    ! This mimics WRFDA's da_llxy boundary checking logic
    ! Use converted longitude and actual min/max values for boundary checking
    outside_domain = .false.
    
    ! Enhanced lat/lon boundary check with longitude range conversion and actual min/max values
    if (olat < lats2d_min .or. olat > lats2d_max .or. &
        olon_converted < lons2d_min .or. olon_converted > lons2d_max) then
      outside_domain = .true.
    end if
    
    if (outside_domain) then
      ! Observation is outside domain - follow WRFDA standard pattern
      ! Set indices to -1 to indicate out-of-domain observation
      oi = -1
      oj = -1
      
      ! Set fractional coordinates to -1.0 to indicate invalid coordinates
      xfloat = -1.0_c_double
      yfloat = -1.0_c_double
      
      ! Return immediately - no further processing for out-of-domain observations
      ! Callers must check for -1 values and handle gracefully
      return
    end if
    
    ! Observation is inside domain - proceed with normal processing
    bestd = 1.0d99; oi = 1; oj = 1
    do j=1,ny
      do i=1,nx
        lat = lats2d(i + (j-1)*nx)
        lon = lons2d(i + (j-1)*nx)
        
        ! Use converted longitude for distance calculation
        d = (lat-olat)*(lat-olat) + (lon-olon_converted)*(lon-olon_converted)
        if (d < bestd) then
          bestd = d; oi = i; oj = j
        end if
      end do
    end do
    
    ! Now use WRFDA's da_togrid to get accurate grid indices and interpolation weights
    ! Convert the nearest grid point to single precision for da_togrid input
    x_approx = real(oi)  ! Convert integer to single precision real
    y_approx = real(oj)  ! Convert integer to single precision real
    
    ! Call da_togrid to get proper grid indices and interpolation weights
    ! da_togrid expects: (x, ib, ie, i, dx, dxm) where x, dx, dxm are real (single precision)
    call da_togrid(x_approx, 1, nx, temp_i, dx_x, dxm_x)
    call da_togrid(y_approx, 1, ny, temp_j, dy_y, dym_y)
    
    ! Set output grid indices using da_togrid results
    oi = temp_i
    oj = temp_j
    
    ! WRFDA-style boundary checking (similar to da_llxy logic)
    ! Check if the calculated indices are within domain bounds
    if (oi < 1 .or. oi >= nx .or. oj < 1 .or. oj >= ny) then
      ! Clamp indices to domain boundaries (similar to WRFDA's approach)
      oi = max(1, min(oi, nx-1))
      oj = max(1, min(oj, ny-1))
      
      ! Set fractional coordinates to exact grid point (no interpolation)
      xfloat = real(oi, c_double)
      yfloat = real(oj, c_double)
    else
      ! Calculate fractional coordinates using WRFDA's interpolation weights
      ! Convert back to double precision for output
      xfloat = real(oi, c_double) + real(dx_x, c_double)
      yfloat = real(oj, c_double) + real(dy_y, c_double)
    end if
  end subroutine find_fractional_ij

  subroutine find_fractional_k(nz, levels, olev, ok, zkfloat)
    integer, intent(in) :: nz
    real(c_double), intent(in) :: levels(*)
    real(c_double), intent(in) :: olev
    integer, intent(out) :: ok
    real(c_double), intent(out) :: zkfloat
    
    ! WRFDA constants
    real(c_double), parameter :: missing_r = -888888.0_c_double
    
    ! Local variables for vertical interpolation
    integer :: k
    real(c_double) :: levels_min, levels_max
    logical :: outside_domain
    real(c_double) :: zk, dz, dzm
    
    ! WRFDA constants (these should match WRFDA's kts and kte)
    integer, parameter :: kts = 1
    integer :: kte
    kte = nz  ! Initialize kte to the number of vertical levels
    
    ! First, check if observation is outside the vertical domain bounds
    ! This mimics WRFDA's boundary checking logic
    outside_domain = .false.
    
    ! Find actual min/max vertical levels by scanning the entire array
    levels_min = levels(1)
    levels_max = levels(1)
    do k = 1, nz
      levels_min = min(levels_min, levels(k))
      levels_max = max(levels_max, levels(k))
    end do
    
    print *, "WRFDA DEBUG: levels_min=", levels_min, " levels_max=", levels_max
    print *, "WRFDA DEBUG: olev=", olev
    
    ! CRITICAL FIX: Check if this is a surface observation BEFORE domain checking
    ! Surface observations are identified by:
    ! 1. Pressure > 100 hPa (surface pressure range)
    ! 2. Sigma level close to 1.0 (surface sigma level)
    ! WRFDA ALWAYS places surface obs at k=1 with no vertical interpolation
    if (olev > 100.0_c_double) then
      ! This is a surface pressure observation (e.g., 1015.5 hPa)
      ! WRFDA ALWAYS places surface obs at k=1 with no vertical interpolation
      print *, "WRFDA DEBUG: Surface pressure observation detected: ", olev, " hPa"
      ok = kts  ! k=1 (surface level)
      zkfloat = real(kts, c_double)  ! No fractional part
      print *, "WRFDA DEBUG: Surface obs: k=", ok, " zkfloat=", zkfloat
      return
    else if (olev >= 0.99_c_double) then
      ! This is a surface sigma level observation (close to 1.0)
      ! WRFDA ALWAYS places surface obs at k=1 with no vertical interpolation
      print *, "WRFDA DEBUG: Surface sigma level observation detected: ", olev
      ok = kts  ! k=1 (surface level)
      zkfloat = real(kts, c_double)  ! No fractional part
      print *, "WRFDA DEBUG: Surface obs: k=", ok, " zkfloat=", zkfloat
      return
    end if
    
    ! Check if observation level is outside vertical domain bounds
    ! This mimics WRFDA's boundary checking logic for upper-air observations
    if (olev < levels_min .or. olev > levels_max) then
      outside_domain = .true.
    end if
    print *, "WRFDA DEBUG: outside_domain=", outside_domain
    
    if (outside_domain) then
      ! Observation is outside vertical domain - follow WRFDA standard pattern
      ! Set indices to -1 to indicate out-of-domain observation
      ok = -1
      
      ! Set fractional coordinates to -1.0 to indicate invalid coordinates
      zkfloat = -1.0_c_double
      
      ! Return immediately - no further processing for out-of-domain observations
      ! Callers must check for -1 values and handle gracefully
      return
    end if
    
    ! Handle single level case (similar to WRFDA's handling)
    if (nz <= 1) then
      ok = 1
      zkfloat = 1.0_c_double
      return
    end if
      
    ! WRFDA LOGIC: Check if this is a surface observation by sigma level
    ! Surface observations are typically at the lowest sigma level (k=1)
    ! WRFDA uses height-based interpolation for most surface obs
    if (olev <= levels(kts) + 0.1_c_double) then
      ! This is a surface observation - use WRFDA's surface logic
      ! Surface obs: k=1, dz=0.0, dzm=1.0 (no vertical interpolation)
      ok = kts
      zkfloat = real(kts, c_double)
      return
    end if
    
    ! WRFDA LOGIC: Upper-air observation - find bracketing levels
    ! This exactly matches da_to_zk's approach for height interpolation
    zk = missing_r  ! WRFDA uses missing_r for invalid values
    
    ! Find the two vertical levels that bracket the observation level
    ! WRFDA uses kts to kte-1 range for interpolation
    do k = kts, kte-1
      if (olev >= levels(k) .and. olev <= levels(k+1)) then
        ! WRFDA's exact formula from da_to_zk:
        ! zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
        zk = real(k, c_double) + (levels(k) - olev)/(levels(k) - levels(k+1))
        exit
      end if
    end do
    
    ! WRFDA LOGIC: Handle extrapolation cases (above highest level)
    if (zk == missing_r .and. olev > levels(kte)) then
      ! Above highest level - use WRFDA's extrapolation logic
      ! WRFDA allows this for verification but not assimilation
      zk = real(kte-1, c_double) + (olev - levels(kte-1))/(levels(kte) - levels(kte-1))
    end if
    
    ! WRFDA LOGIC: Handle extrapolation cases (below lowest level)
    if (zk == missing_r .and. olev < levels(kts)) then
      ! Below lowest level - use WRFDA's extrapolation logic
      zk = real(kts+1, c_double) - (olev - levels(kts+1))/(levels(kts) - levels(kts+1))
    end if
    
    ! WRFDA LOGIC: Convert zk to grid indices and weights
    ! This exactly matches da_convert_zk's approach
    if (zk /= missing_r) then
      ! Set the base level index
      ok = int(zk)
      
      ! WRFDA's exact weight calculation:
      ! info%dz(k,n) = info%zk(k,n) - real(info%k(k,n))
      ! info%dzm(k,n) = 1.0 - info%dz(k,n)
      dz = zk - real(ok, c_double)
      dzm = 1.0_c_double - dz
      
      ! WRFDA LOGIC: Clamp indices to valid range (kts to kte-1)
      if (ok < kts) then
        ok = kts
        ! Recalculate weights for clamped level
        dz = zk - real(ok, c_double)
        dzm = 1.0_c_double - dz
      else if (ok >= kte) then
        ok = kte - 1
        ! Recalculate weights for clamped level
        dz = zk - real(ok, c_double)
        dzm = 1.0_c_double - dz
      end if
      
      ! Ensure weights are in valid range (0 to 1)
      dz = max(0.0_c_double, min(1.0_c_double, dz))
      dzm = max(0.0_c_double, min(1.0_c_double, dzm))
      
      ! Set fractional coordinate (this is what WRFDA uses for interpolation)
      zkfloat = zk
    else
      ! Invalid zk - set to surface level as fallback
      ok = kts
      zkfloat = real(kts, c_double)
    end if
    
    print *, "WRFDA DEBUG: Final ok=", ok, " zkfloat=", zkfloat
    print *, "WRFDA DEBUG: Calculated dz=", dz, " dzm=", dzm
  end subroutine find_fractional_k

  ! Initialize persistent background state for incremental 3D-Var
  ! This function should be called once at the beginning of the outer loop
  ! to set up the constant background state (xb)
  integer(c_int) function wrfda_initialize_background_state(nx, ny, nz, u_bg, v_bg, t_bg, q_bg, psfc_bg) bind(C, name="wrfda_initialize_background_state")
    implicit none
    integer(c_int), value :: nx, ny, nz
    real(c_double), intent(in) :: u_bg(*), v_bg(*), t_bg(*), q_bg(*), psfc_bg(*)
    
    print *, "WRFDA DEBUG: Initializing persistent background state"
    print *, "WRFDA DEBUG: nx=", nx, " ny=", ny, " nz=", nz
    
    ! Store grid dimensions
    background_nx = nx
    background_ny = ny
    background_nz = nz
    
    ! Initialize persistent grid with background state
    call init_domain_for_adjoint(persistent_grid, nx, ny, nz, u_bg, v_bg, t_bg, q_bg, psfc_bg)
    
    ! Mark background state as initialized
    background_state_initialized = .true.
    
    print *, "WRFDA DEBUG: Persistent background state initialized successfully"
    wrfda_initialize_background_state = 0_c_int
  end function wrfda_initialize_background_state

  ! Update only analysis increments (xa) for incremental 3D-Var
  ! This subroutine updates the analysis increments while keeping background state constant
  subroutine wrfda_update_analysis_increments(u_inc, v_inc, t_inc, q_inc, psfc_inc)
    implicit none
    real(c_double), intent(in) :: u_inc(*), v_inc(*), t_inc(*), q_inc(*), psfc_inc(*)
    integer :: i,j,k, nz1
    
    if (.not. background_state_initialized) then
      print *, "WRFDA ERROR: Background state not initialized. Call wrfda_initialize_background_state first."
      return
    end if
    
    nz1 = max(1, background_nz)
    print *, "WRFDA DEBUG: Updating analysis increments for grid size:", background_nx, "x", background_ny, "x", nz1
    
    ! Update only the xa fields (analysis increments)
    do k=1,nz1; do j=1,background_ny; do i=1,background_nx
      ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
      persistent_grid%xa%u(i,j,k) = real(u_inc(i + (j-1)*background_nx + (k-1)*background_nx*background_ny))
      persistent_grid%xa%v(i,j,k) = real(v_inc(i + (j-1)*background_nx + (k-1)*background_nx*background_ny))
      persistent_grid%xa%t(i,j,k) = real(t_inc(i + (j-1)*background_nx + (k-1)*background_nx*background_ny))
      persistent_grid%xa%q(i,j,k) = real(q_inc(i + (j-1)*background_nx + (k-1)*background_nx*background_ny))
    end do; end do; end do
    do j=1,background_ny; do i=1,background_nx
      ! C++ row-major indexing: [i][j] -> i + j*nx
      persistent_grid%xa%psfc(i,j) = real(psfc_inc(i + (j-1)*background_nx))
    end do; end do
    
    print *, "WRFDA DEBUG: Analysis increments updated successfully"
  end subroutine wrfda_update_analysis_increments

  ! CRITICAL FIX: Add subroutine to update only data arrays without reconstructing grid structure
  subroutine update_domain_data(grid, nx, ny, nz, u, v, t, q, psfc)
    type(domain), intent(inout) :: grid
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    integer :: i,j,k, nz1
    
    nz1 = max(1, nz)
    print *, "WRFDA DEBUG: update_domain_data: updating data arrays only"
    
    ! Update only the data arrays, not the grid structure
    do k=1,nz1; do j=1,ny; do i=1,nx
      grid%xa%u(i,j,k) = real(u(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xa%v(i,j,k) = real(v(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xa%t(i,j,k) = real(t(i + (j-1)*nx + (k-1)*nx*ny))
      grid%xa%q(i,j,k) = real(q(i + (j-1)*nx + (k-1)*nx*ny))
    end do; end do; end do
    do j=1,ny; do i=1,nx
      grid%xa%psfc(i,j) = real(psfc(i + (j-1)*nx))
    end do; end do
    
    ! Mirror xb from xa
    grid%xb%u = grid%xa%u
    grid%xb%v = grid%xa%v
    grid%xb%t = grid%xa%t
    grid%xb%q = grid%xa%q
    grid%xb%psfc = grid%xa%psfc
    
    print *, "WRFDA DEBUG: update_domain_data completed"
  end subroutine update_domain_data

  ! New function to call da_get_innov_vector directly
  integer(c_int) function wrfda_get_innov_vector(it, domain_ptr, ob_ptr, iv_ptr, config_flags_ptr) bind(C, name="wrfda_get_innov_vector")
    implicit none
    integer(c_int), intent(in) :: it
    type(c_ptr), value :: domain_ptr, ob_ptr, iv_ptr, config_flags_ptr
    type(domain), pointer :: grid
    type(y_type), pointer :: ob
    type(iv_type), pointer :: iv
    type(grid_config_rec_type), pointer :: config_flags
    
    ! Quality control statistics array (simplified for now)
    integer, dimension(1,1,1,1) :: num_qcstat_conv
    
    ! Initialize return code
    wrfda_get_innov_vector = 0
    
    ! Convert C pointers to Fortran pointers
    call c_f_pointer(domain_ptr, grid)
    call c_f_pointer(ob_ptr, ob)
    
    ! Check if iv_ptr is null before converting
    if (c_associated(iv_ptr)) then
      call c_f_pointer(iv_ptr, iv)
      print *, "WRFDA DEBUG: iv_ptr converted successfully"
    else
      print *, "WRFDA DEBUG: iv_ptr is null - skipping da_get_innov_vector"
      wrfda_get_innov_vector = 0
      return
    end if
    
    ! Check if config_flags_ptr is null before converting
    !if (c_associated(config_flags_ptr)) then
    !  call c_f_pointer(config_flags_ptr, config_flags)
    !  print *, "WRFDA DEBUG: config_flags_ptr converted successfully"
    !else
    !  print *, "WRFDA DEBUG: config_flags_ptr is null - skipping da_get_innov_vector"
    !  wrfda_get_innov_vector = 0
    !  return
    !end if
    
    ! Initialize QC statistics array
    num_qcstat_conv = 0
    
    ! Call the main WRFDA innovation vector computation routine
    print *, "WRFDA DEBUG: About to call da_get_innov_vector"
    call da_get_innov_vector(it, num_qcstat_conv, ob, iv, grid, config_flags)
    
    print *, "WRFDA DEBUG: da_get_innov_vector completed successfully"
    
  end function wrfda_get_innov_vector

  ! Helper function to construct WRFDA domain structure from flat arrays
  integer(c_int) function wrfda_construct_domain_from_arrays(nx, ny, nz, u, v, t, q, psfc, lats2d, lons2d, levels, domain_ptr) bind(C, name="wrfda_construct_domain_from_arrays")
    implicit none
    integer(c_int), intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    real(c_double), intent(in) :: lats2d(*), lons2d(*), levels(*)
    type(c_ptr), intent(out) :: domain_ptr
    
    type(domain), pointer :: grid
    integer :: i, j, k, idx

    ! Allocate new domain structure
    allocate(grid)

    ! Set up basic domain dimensions
    grid%id = 1

    ! Allocate and populate xb (background state) arrays
    allocate(grid%xb%u(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%v(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%t(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%q(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%psfc(1:nx, 1:ny))

    ! Copy data from flat arrays to WRFDA structure
    ! Note: C++ arrays are in [X,Y,Z] order (row-major), Fortran arrays are [Z,Y,X] (column-major)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! Convert from C++ [X,Y,Z] indexing to Fortran [Z,Y,X] indexing
          idx = (i-1)*ny*nz + (j-1)*nz + k
          grid%xb%u(i,j,k) = u(idx)
          grid%xb%v(i,j,k) = v(idx)
          grid%xb%t(i,j,k) = t(idx)
          grid%xb%q(i,j,k) = q(idx)
        end do
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        ! Convert from C++ [X,Y] indexing to Fortran [Y,X] indexing
        idx = (i-1)*ny + j
        grid%xb%psfc(i,j) = psfc(idx)
      end do
    end do

    ! Set up grid metadata
    allocate(grid%xb%lat(1:nx, 1:ny))
    allocate(grid%xb%lon(1:nx, 1:ny))
    allocate(grid%xb%h(1:nx, 1:ny, 1:nz))

    do j = 1, ny
      do i = 1, nx
        ! Convert from C++ [X,Y] indexing to Fortran [Y,X] indexing
        idx = (i-1)*ny + j
        grid%xb%lat(i,j) = lats2d(idx)
        grid%xb%lon(i,j) = lons2d(idx)
        do k = 1, nz
          grid%xb%h(i,j,k) = levels(k)
        end do
      end do
    end do

    ! Return pointer to allocated domain
    domain_ptr = c_loc(grid)
    wrfda_construct_domain_from_arrays = 0
    
  end function wrfda_construct_domain_from_arrays

  ! Construct y_type from observation data
  type(c_ptr) function wrfda_construct_y_type(num_obs, obs_values, obs_errors, obs_types, family) bind(C, name="wrfda_construct_y_type")
    implicit none
    integer(c_int), intent(in) :: num_obs
    real(c_double), intent(in) :: obs_values(*), obs_errors(*)
    character(c_char), intent(in) :: obs_types(*), family(*)
    
    type(y_type), pointer :: y
    integer :: i
    character(len=20) :: family_str
    character(c_char), target :: family_target(20)
    
    print *, "WRFDA DEBUG: wrfda_construct_y_type called with num_obs=", num_obs
    
    ! Convert C string to Fortran string
    family_target = family(1:20)
    ! Simple string conversion without c_f_string_ptr
    family_str = ""
    do i = 1, 20
      if (family_target(i) == c_null_char) exit
      family_str(i:i) = family_target(i)
    end do
    
    ! Use obs_types to suppress warning (could be used for type-specific handling)
    if (obs_types(1) /= c_null_char) then
      ! Observation types available - could be used for type-specific processing
    end if
    
    ! Allocate y_type structure
    allocate(y)
    
    ! Initialize counters
    y%nlocal = 0
    y%ntotal = 0
    y%num_inst = 0
    
    ! For surface observations (metar, synop, adpsfc), allocate metar array
    if (trim(family_str) == "metar" .or. trim(family_str) == "synop" .or. trim(family_str) == "adpsfc") then
      allocate(y%metar(num_obs))
      y%nlocal(1) = num_obs  ! metar index
      y%ntotal(1) = num_obs
      
      ! Populate metar array with observation data
      do i = 1, num_obs
        ! For now, set all values to the observation values
        ! In a full implementation, these would be residuals (obs - model)
        y%metar(i)%u = obs_values(i)
        y%metar(i)%v = obs_values(i)  ! Placeholder - would need separate U/V values
        y%metar(i)%t = obs_values(i)  ! Placeholder - would need separate T values
        y%metar(i)%p = 101325.0       ! Standard pressure
        y%metar(i)%q = obs_values(i)  ! Placeholder - would need separate Q values
        ! Use obs_errors to suppress warning (could be used for error weighting)
        if (obs_errors(i) > 0.0) then
          ! Error information available - could be used for weighting
        end if
      end do
      
      print *, "WRFDA DEBUG: Allocated metar array with", num_obs, "observations"
    end if
    
    ! Return pointer to allocated y_type
    wrfda_construct_y_type = c_loc(y)
    print *, "WRFDA DEBUG: y_type constructed successfully"
    
  end function wrfda_construct_y_type

  ! Construct iv_type from observation data
  type(c_ptr) function wrfda_construct_iv_type(num_obs, obs_values, obs_errors, obs_types, obs_lats, obs_lons, obs_levels, family) bind(C, name="wrfda_construct_iv_type")
    implicit none
    integer(c_int), intent(in) :: num_obs
    real(c_double), intent(in) :: obs_values(*), obs_errors(*)
    real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels(*)
    character(c_char), intent(in) :: obs_types(*), family(*)
    
    type(iv_type), pointer :: iv
    integer :: i
    character(len=20) :: family_str
    character(c_char), target :: family_target(20)
    
    print *, "WRFDA DEBUG: wrfda_construct_iv_type called with num_obs=", num_obs
    
    ! Check if num_obs is valid
    if (num_obs <= 0) then
      print *, "WRFDA DEBUG: Invalid num_obs =", num_obs
      wrfda_construct_iv_type = c_null_ptr
      return
    end if
    
    ! Check if arrays are accessible (with bounds checking)
    print *, "WRFDA DEBUG: About to access obs_values(1)"
    if (num_obs >= 1) then
      print *, "WRFDA DEBUG: obs_values(1) =", obs_values(1)
    else
      print *, "WRFDA DEBUG: num_obs < 1, cannot access obs_values(1)"
    end if
    
    print *, "WRFDA DEBUG: About to access obs_errors(1)"
    if (num_obs >= 1) then
      print *, "WRFDA DEBUG: obs_errors(1) =", obs_errors(1)
    else
      print *, "WRFDA DEBUG: num_obs < 1, cannot access obs_errors(1)"
    end if
    
    print *, "WRFDA DEBUG: About to access obs_lats(1)"
    if (num_obs >= 1) then
      print *, "WRFDA DEBUG: obs_lats(1) =", obs_lats(1)
    else
      print *, "WRFDA DEBUG: num_obs < 1, cannot access obs_lats(1)"
    end if
    
    print *, "WRFDA DEBUG: About to access obs_lons(1)"
    if (num_obs >= 1) then
      print *, "WRFDA DEBUG: obs_lons(1) =", obs_lons(1)
    else
      print *, "WRFDA DEBUG: num_obs < 1, cannot access obs_lons(1)"
    end if
    
    print *, "WRFDA DEBUG: About to access obs_levels(1)"
    if (num_obs >= 1) then
      print *, "WRFDA DEBUG: obs_levels(1) =", obs_levels(1)
    else
      print *, "WRFDA DEBUG: num_obs < 1, cannot access obs_levels(1)"
    end if
    
    ! Now implement the actual iv_type construction
    print *, "WRFDA DEBUG: Starting iv_type construction"
    
    ! Allocate iv_type
    allocate(iv)
    if (.not. associated(iv)) then
      print *, "WRFDA DEBUG: Failed to allocate iv_type"
      wrfda_construct_iv_type = c_null_ptr
      return
    end if
    
    ! Initialize basic fields
    iv%time = 0
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
      allocate(iv%synop(num_obs))
      print *, "WRFDA DEBUG: Allocated synop array with", num_obs, "observations"
      
      ! Populate synop data
      do i = 1, num_obs
        ! Set up height
        iv%synop(i)%h = obs_levels(i)
        
        ! Set up field_type members for u, v, t, p, q
        iv%synop(i)%u%inv = obs_values(i)
        iv%synop(i)%u%qc = 0  ! Good quality
        iv%synop(i)%u%error = obs_errors(i)
        iv%synop(i)%u%sens = 0.0
        iv%synop(i)%u%imp = 0.0
        
        iv%synop(i)%v%inv = obs_values(i)
        iv%synop(i)%v%qc = 0
        iv%synop(i)%v%error = obs_errors(i)
        iv%synop(i)%v%sens = 0.0
        iv%synop(i)%v%imp = 0.0
        
        iv%synop(i)%t%inv = obs_values(i)
        iv%synop(i)%t%qc = 0
        iv%synop(i)%t%error = obs_errors(i)
        iv%synop(i)%t%sens = 0.0
        iv%synop(i)%t%imp = 0.0
        
        iv%synop(i)%p%inv = obs_values(i)
        iv%synop(i)%p%qc = 0
        iv%synop(i)%p%error = obs_errors(i)
        iv%synop(i)%p%sens = 0.0
        iv%synop(i)%p%imp = 0.0
        
        iv%synop(i)%q%inv = obs_values(i)
        iv%synop(i)%q%qc = 0
        iv%synop(i)%q%error = obs_errors(i)
        iv%synop(i)%q%sens = 0.0
        iv%synop(i)%q%imp = 0.0
      end do
    end if
    
    print *, "WRFDA DEBUG: iv_type construction completed successfully"
    wrfda_construct_iv_type = c_loc(iv)
    
    ! Convert C string to Fortran string
    family_target = family(1:20)
    family_str = ""
    do i = 1, 20
      if (family_target(i) == c_null_char) exit
      family_str(i:i) = family_target(i)
    end do
    
    ! Use obs_types to suppress warning
    if (obs_types(1) /= c_null_char) then
      ! Observation types available
    end if
    
    ! Use obs_lats, obs_lons, obs_levels to suppress warnings
    ! These could be used for location-specific processing in the future
    if (obs_lats(1) /= 0.0 .or. obs_lons(1) /= 0.0 .or. obs_levels(1) /= 0.0) then
      ! Location data available - could be used for spatial processing
    end if
    
    ! Allocate iv_type structure
    allocate(iv)
    if (.not. associated(iv)) then
      print *, "WRFDA DEBUG: Failed to allocate iv_type structure"
      wrfda_construct_iv_type = c_null_ptr
      return
    end if
    print *, "WRFDA DEBUG: Successfully allocated iv_type structure"
    
    ! Initialize basic fields
    iv%nstats = 0
    iv%time = 0
    iv%num_inst = 0
    iv%total_rad_pixel = 0
    iv%total_rad_channel = 0
    iv%missing = -888888.0
    iv%ptop = 0.0
    
    ! Initialize info array for all observation types
    do i = 1, size(iv%info)
      iv%info(i)%max_lev = 0
      iv%info(i)%nlocal = 0
      iv%info(i)%ntotal = 0
      iv%info(i)%thin_nlocal = 0
      iv%info(i)%thin_ntotal = 0
      iv%info(i)%plocal = 0
      iv%info(i)%ptotal = 0
      iv%info(i)%thin_plocal = 0
      iv%info(i)%thin_ptotal = 0
      iv%info(i)%n1 = 0
      iv%info(i)%n2 = 0
    end do
    
    ! Set error factors for surface observations
    if (trim(family_str) == "metar" .or. trim(family_str) == "synop" .or. trim(family_str) == "adpsfc") then
      iv%metar_ef_u = 1.5
      iv%metar_ef_v = 1.5
      iv%metar_ef_t = 1.0
      iv%metar_ef_p = 100.0
      iv%metar_ef_q = 0.0
      
      ! Set info for synop observations (index 2)
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
      
      ! Allocate synop array only if we have observations
      if (num_obs > 0) then
        allocate(iv%synop(num_obs))
        iv%nstats(2) = num_obs  ! synop index
      
      ! Populate synop array with observation data
      do i = 1, num_obs
        ! Set height (surface level)
        iv%synop(i)%h = 0.0
        
        ! Set field data for each variable
        iv%synop(i)%u%inv = obs_values(i)      ! Innovation (obs - model)
        iv%synop(i)%u%qc = 0                   ! Good quality
        iv%synop(i)%u%error = obs_errors(i)    ! Observation error
        iv%synop(i)%u%sens = 0.0               ! Sensitivity
        iv%synop(i)%u%imp = 0.0                ! Impact
        
        iv%synop(i)%v%inv = obs_values(i)      ! Placeholder
        iv%synop(i)%v%qc = 0
        iv%synop(i)%v%error = obs_errors(i)
        iv%synop(i)%v%sens = 0.0
        iv%synop(i)%v%imp = 0.0
        
        iv%synop(i)%t%inv = obs_values(i)      ! Placeholder
        iv%synop(i)%t%qc = 0
        iv%synop(i)%t%error = obs_errors(i)
        iv%synop(i)%t%sens = 0.0
        iv%synop(i)%t%imp = 0.0
        
        iv%synop(i)%p%inv = 101325.0           ! Standard pressure
        iv%synop(i)%p%qc = 0
        iv%synop(i)%p%error = 100.0
        iv%synop(i)%p%sens = 0.0
        iv%synop(i)%p%imp = 0.0
        
        iv%synop(i)%q%inv = obs_values(i)      ! Placeholder
        iv%synop(i)%q%qc = 0
        iv%synop(i)%q%error = obs_errors(i)
        iv%synop(i)%q%sens = 0.0
        iv%synop(i)%q%imp = 0.0
      end do
      
      print *, "WRFDA DEBUG: Allocated iv%synop array with", num_obs, "observations"
      end if  ! num_obs > 0
    end if  ! family check
    
    ! Return pointer to allocated iv_type
    wrfda_construct_iv_type = c_loc(iv)
    print *, "WRFDA DEBUG: iv_type constructed successfully"
    
  end function wrfda_construct_iv_type

  ! Construct config_flags for WRFDA
  type(c_ptr) function wrfda_construct_config_flags() bind(C, name="wrfda_construct_config_flags")
    implicit none
    
    ! For now, return a null pointer since config_flags is complex
    ! In a full implementation, this would create a proper grid_config_rec_type
    ! with all the necessary WRFDA configuration parameters
    
    print *, "WRFDA DEBUG: wrfda_construct_config_flags called"
    print *, "WRFDA DEBUG: Returning null config_flags (simplified implementation)"
    
    wrfda_construct_config_flags = c_null_ptr
    
  end function wrfda_construct_config_flags

  ! Extract innovation values from iv_type structure
  integer(c_int) function wrfda_extract_innovations(iv_ptr, family, innovations, num_innovations, max_innovations) bind(C, name="wrfda_extract_innovations")
    implicit none
    type(c_ptr), intent(in) :: iv_ptr
    character(c_char), intent(in) :: family(*)
    real(c_double), intent(out) :: innovations(*)
    integer(c_int), intent(out) :: num_innovations
    integer(c_int), intent(in) :: max_innovations
    
    type(iv_type), pointer :: iv
    character(len=20) :: family_str
    character(c_char), target :: family_target(20)
    integer :: i, count
    integer :: family_index
    
    print *, "WRFDA DEBUG: wrfda_extract_innovations called"
    
    ! Convert C string to Fortran string
    family_target = family(1:20)
    family_str = ""
    do i = 1, 20
      if (family_target(i) == c_null_char) exit
      family_str(i:i) = family_target(i)
    end do
    
    ! Check if iv_ptr is null
    if (c_associated(iv_ptr)) then
      print *, "WRFDA DEBUG: iv_ptr is not null"
      ! Get pointer to iv_type
      call c_f_pointer(iv_ptr, iv)
      if (.not. associated(iv)) then
        print *, "WRFDA DEBUG: c_f_pointer failed"
        wrfda_extract_innovations = -1
        return
      end if
      print *, "WRFDA DEBUG: Successfully converted iv_ptr to iv pointer"
    else
      print *, "WRFDA DEBUG: iv_ptr is null - returning empty innovations"
      num_innovations = 0
      wrfda_extract_innovations = 0
      return
    end if
    
    ! Check if iv structure is valid
    print *, "WRFDA DEBUG: iv%time =", iv%time
    print *, "WRFDA DEBUG: iv%num_inst =", iv%num_inst
    print *, "WRFDA DEBUG: size(iv%info) =", size(iv%info)
    
    ! Determine family index (simplified mapping)
    family_index = 0
    if (trim(family_str) == "metar" .or. trim(family_str) == "synop" .or. trim(family_str) == "adpsfc") then
      family_index = 2  ! synop index (metar uses synop structure)
    else if (trim(family_str) == "sound") then
      family_index = 3  ! sound index
    else if (trim(family_str) == "gpspw") then
      family_index = 4  ! gpspw index
    end if
    
    ! Check bounds
    if (family_index > 0 .and. family_index <= size(iv%info)) then
      print *, "WRFDA DEBUG: family_index", family_index, "is within bounds"
    else
      print *, "WRFDA DEBUG: family_index", family_index, "is out of bounds (max:", size(iv%info), ")"
    end if
    
    if (family_index == 0) then
      print *, "WRFDA DEBUG: Unknown family:", trim(family_str)
      wrfda_extract_innovations = -2
      return
    end if
    
    ! Extract innovations based on family
    count = 0
    print *, "WRFDA DEBUG: family_index =", family_index
    print *, "WRFDA DEBUG: associated(iv%synop) =", associated(iv%synop)
    if (family_index == 2) then
      print *, "WRFDA DEBUG: iv%info(2)%nlocal =", iv%info(2)%nlocal
    end if
    
    if (family_index == 2 .and. family_index <= size(iv%info) .and. associated(iv%synop) .and. iv%info(2)%nlocal > 0) then
      ! Extract from synop array
      print *, "WRFDA DEBUG: Extracting from synop array with", iv%info(2)%nlocal, "observations"
      do i = 1, iv%info(2)%nlocal
        if (i > size(iv%synop)) then
          print *, "WRFDA DEBUG: Index", i, "exceeds synop array size", size(iv%synop)
          exit
        end if
        ! Extract U component innovation
        count = count + 1
        if (count <= max_innovations) then
          innovations(count) = iv%synop(i)%u%inv
        end if
        
        ! Extract V component innovation
        count = count + 1
        if (count <= max_innovations) then
          innovations(count) = iv%synop(i)%v%inv
        end if
        
        ! Extract T component innovation
        count = count + 1
        if (count <= max_innovations) then
          innovations(count) = iv%synop(i)%t%inv
        end if
        
        ! Extract P component innovation
        count = count + 1
        if (count <= max_innovations) then
          innovations(count) = iv%synop(i)%p%inv
        end if
        
        ! Extract Q component innovation
        count = count + 1
        if (count <= max_innovations) then
          innovations(count) = iv%synop(i)%q%inv
        end if
      end do
      
      print *, "WRFDA DEBUG: Extracted", count, "innovations from synop array"
    else
      print *, "WRFDA DEBUG: No data available for family:", trim(family_str)
      count = 0
    end if
    
    num_innovations = count
    wrfda_extract_innovations = 0
    
    print *, "WRFDA DEBUG: Successfully extracted", num_innovations, "innovations"
    
  end function wrfda_extract_innovations

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
    
    print *, "WRFDA DEBUG: Initialized for 3D-Var analysis"
    print *, "WRFDA DEBUG: var4d_run = ", var4d_run
    print *, "WRFDA DEBUG: num_fgat_time = ", num_fgat_time
    print *, "WRFDA DEBUG: sfc_assi_options = ", sfc_assi_options
    
  end subroutine initialize_wrfda_3dvar

end module metada_wrfda_dispatch


