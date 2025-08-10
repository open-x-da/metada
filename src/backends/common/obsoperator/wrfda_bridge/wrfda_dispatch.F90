! Dispatch C-bindings that translate flat arrays into WRFDA types
! and call the included transform routines.

module metada_wrfda_dispatch
  use iso_c_binding
  use module_domain,        only: domain, x_type
  use da_define_structures, only: iv_type, y_type, da_allocate_y
  use da_control, only: metar, synop, ships, buoy, airep, pilot, sound, sonde_sfc, &
                        sfc_assi_options, sfc_assi_options_1, trace_use_dull
  use da_metar,  only: da_transform_xtoy_metar,  da_transform_xtoy_metar_adj
  use da_synop,  only: da_transform_xtoy_synop,  da_transform_xtoy_synop_adj
  use da_buoy,   only: da_transform_xtoy_buoy,   da_transform_xtoy_buoy_adj
  use da_ships,  only: da_transform_xtoy_ships,  da_transform_xtoy_ships_adj
  use da_airep,  only: da_transform_xtoy_airep,  da_transform_xtoy_airep_adj
  use da_pilot,  only: da_transform_xtoy_pilot,  da_transform_xtoy_pilot_adj
  use da_sound,  only: da_transform_xtoy_sound,  da_transform_xtoy_sound_adj, &
                       da_transform_xtoy_sonde_sfc, da_transform_xtoy_sonde_sfc_adj
  implicit none
contains

  ! For brevity, we'll implement a simplified bridge that maps only a uniform field
  ! to observations via WRFDA calls is non-trivial without full grid/iv construction
  ! from our C arrays. Here we expose the C symbols but currently return error until
  ! a proper mapping is implemented.

  integer(c_int) function wrfda_xtoy_apply(operator_family, wrfda_root, state_values, nx, ny, nz, obs_lats, obs_lons, obs_levels, num_obs, out_y) bind(C, name="wrfda_xtoy_apply")
    implicit none
    character(c_char), intent(in) :: operator_family(*)
    character(c_char), intent(in) :: wrfda_root(*)
    real(c_double), intent(in) :: state_values(*)
    integer(c_int), value :: nx, ny, nz
    real(c_double), intent(in) :: obs_lats(*)
    real(c_double), intent(in) :: obs_lons(*)
    real(c_double), intent(in) :: obs_levels(*)
    integer(c_int), value :: num_obs
    real(c_double), intent(out) :: out_y(*)
    integer :: n
    integer :: inx, iny, inz
    real(c_double) :: x, y, z
    integer :: i0, i1, j0, j1, k0, k1
    real(c_double) :: wx0, wx1, wy0, wy1, wz0, wz1
    real(c_double) :: val
    integer :: idx000, idx100, idx010, idx110, idx001, idx101, idx011, idx111

    inx = max(1, nx); iny = max(1, ny); inz = max(1, nz)

    do n = 1, num_obs
      ! Interpret obs_lats/obs_lons/obs_levels as fractional grid indices (1-based)
      x = max(1.0_c_double, min(real(inx, c_double), obs_lats(n)))
      y = max(1.0_c_double, min(real(iny, c_double), obs_lons(n)))
      z = max(1.0_c_double, min(real(max(inz,1), c_double), merge(obs_levels(n), 1.0_c_double, inz>1)))

      i0 = int(floor(x)); i1 = min(i0 + 1, inx)
      j0 = int(floor(y)); j1 = min(j0 + 1, iny)
      k0 = int(floor(z)); k1 = min(k0 + 1, max(inz,1))

      wx1 = x - real(i0, c_double); wx0 = 1.0_c_double - wx1
      wy1 = y - real(j0, c_double); wy0 = 1.0_c_double - wy1
      if (inz > 1) then
        wz1 = z - real(k0, c_double); wz0 = 1.0_c_double - wz1
      else
        wz0 = 1.0_c_double; wz1 = 0.0_c_double
        k0 = 1; k1 = 1
      end if

      ! Convert to linear indices (Fortran 1-based, column-major: i fastest)
      idx000 = i0 + (j0-1)*inx + (k0-1)*inx*iny
      idx100 = i1 + (j0-1)*inx + (k0-1)*inx*iny
      idx010 = i0 + (j1-1)*inx + (k0-1)*inx*iny
      idx110 = i1 + (j1-1)*inx + (k0-1)*inx*iny
      idx001 = i0 + (j0-1)*inx + (k1-1)*inx*iny
      idx101 = i1 + (j0-1)*inx + (k1-1)*inx*iny
      idx011 = i0 + (j1-1)*inx + (k1-1)*inx*iny
      idx111 = i1 + (j1-1)*inx + (k1-1)*inx*iny

      val = &
        wx0*wy0*wz0*state_values(idx000) + &
        wx1*wy0*wz0*state_values(idx100) + &
        wx0*wy1*wz0*state_values(idx010) + &
        wx1*wy1*wz0*state_values(idx110) + &
        wx0*wy0*wz1*state_values(idx001) + &
        wx1*wy0*wz1*state_values(idx101) + &
        wx0*wy1*wz1*state_values(idx011) + &
        wx1*wy1*wz1*state_values(idx111)

      out_y(n) = val
    end do
    wrfda_xtoy_apply = 0_c_int
  end function wrfda_xtoy_apply

  integer(c_int) function wrfda_xtoy_adjoint(operator_family, wrfda_root, delta_y, obs_lats, obs_lons, obs_levels, num_obs, inout_state_values, nx, ny, nz) bind(C, name="wrfda_xtoy_adjoint")
    implicit none
    character(c_char), intent(in) :: operator_family(*)
    character(c_char), intent(in) :: wrfda_root(*)
    real(c_double), intent(in) :: delta_y(*)
    real(c_double), intent(in) :: obs_lats(*)
    real(c_double), intent(in) :: obs_lons(*)
    real(c_double), intent(in) :: obs_levels(*)
    integer(c_int), value :: num_obs
    real(c_double), intent(inout) :: inout_state_values(*)
    integer(c_int), value :: nx, ny, nz
    integer :: n
    integer :: inx, iny, inz
    real(c_double) :: x, y, z
    integer :: i0, i1, j0, j1, k0, k1
    real(c_double) :: wx1, wx0, wy1, wy0, wz1, wz0
    real(c_double) :: w
    integer :: idx000, idx100, idx010, idx110, idx001, idx101, idx011, idx111

    inx = max(1, nx); iny = max(1, ny); inz = max(1, nz)

    do n = 1, num_obs
      x = max(1.0_c_double, min(real(inx, c_double), obs_lats(n)))
      y = max(1.0_c_double, min(real(iny, c_double), obs_lons(n)))
      z = max(1.0_c_double, min(real(max(inz,1), c_double), merge(obs_levels(n), 1.0_c_double, inz>1)))

      i0 = int(floor(x)); i1 = min(i0 + 1, inx)
      j0 = int(floor(y)); j1 = min(j0 + 1, iny)
      k0 = int(floor(z)); k1 = min(k0 + 1, max(inz,1))

      wx1 = x - real(i0, c_double); wx0 = 1.0_c_double - wx1
      wy1 = y - real(j0, c_double); wy0 = 1.0_c_double - wy1
      if (inz > 1) then
        wz1 = z - real(k0, c_double); wz0 = 1.0_c_double - wz1
      else
        wz0 = 1.0_c_double; wz1 = 0.0_c_double
        k0 = 1; k1 = 1
      end if

      idx000 = i0 + (j0-1)*inx + (k0-1)*inx*iny
      idx100 = i1 + (j0-1)*inx + (k0-1)*inx*iny
      idx010 = i0 + (j1-1)*inx + (k0-1)*inx*iny
      idx110 = i1 + (j1-1)*inx + (k0-1)*inx*iny
      idx001 = i0 + (j0-1)*inx + (k1-1)*inx*iny
      idx101 = i1 + (j0-1)*inx + (k1-1)*inx*iny
      idx011 = i0 + (j1-1)*inx + (k1-1)*inx*iny
      idx111 = i1 + (j1-1)*inx + (k1-1)*inx*iny

      ! Distribute adjoint increment to neighbors
      w = wx0*wy0*wz0; inout_state_values(idx000) = inout_state_values(idx000) + w*delta_y(n)
      w = wx1*wy0*wz0; inout_state_values(idx100) = inout_state_values(idx100) + w*delta_y(n)
      w = wx0*wy1*wz0; inout_state_values(idx010) = inout_state_values(idx010) + w*delta_y(n)
      w = wx1*wy1*wz0; inout_state_values(idx110) = inout_state_values(idx110) + w*delta_y(n)
      w = wx0*wy0*wz1; inout_state_values(idx001) = inout_state_values(idx001) + w*delta_y(n)
      w = wx1*wy0*wz1; inout_state_values(idx101) = inout_state_values(idx101) + w*delta_y(n)
      w = wx0*wy1*wz1; inout_state_values(idx011) = inout_state_values(idx011) + w*delta_y(n)
      w = wx1*wy1*wz1; inout_state_values(idx111) = inout_state_values(idx111) + w*delta_y(n)
    end do

    wrfda_xtoy_adjoint = 0_c_int
  end function wrfda_xtoy_adjoint

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

  ! Real WRFDA call via grid arrays: build minimal domain/iv/y and dispatch
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

    type(domain), target :: grid
    type(iv_type), target :: iv
    type(y_type), target :: y
    integer :: n
    character(len=256) :: fambuf
    integer :: famlen
    character(len=:), allocatable :: op_str, fam_str, var_str
    integer :: fam_id
    character(len=1) :: var_code

    ! Minimal init (high-level: we assume single tile, column-major layout)
    call init_domain_from_arrays(grid, nx, ny, nz, u, v, t, q, psfc)
    sfc_assi_options = sfc_assi_options_1

    fambuf = '' ; famlen = 0
    do n = 1, len(fambuf)
      if (operator_family(n) == c_null_char) exit
      famlen = famlen + 1
      fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
    end do
    op_str = trim(fambuf(1:max(famlen,1)))
    call parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
    call init_iv_from_obs(fam_id, iv, nx, ny, nz, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels)
    call da_allocate_y(iv, y)

    select case (trim(fam_str))
    case ('metar'); call da_transform_xtoy_metar(grid, iv, y)
    case ('synop'); call da_transform_xtoy_synop(grid, iv, y)
    case ('buoy');  call da_transform_xtoy_buoy(grid, iv, y)
    case ('ships'); call da_transform_xtoy_ships(grid, iv, y)
    case ('airep'); call da_transform_xtoy_airep(grid, iv, y)
    case ('pilot'); call da_transform_xtoy_pilot(grid, iv, y)
    case ('sound'); call da_transform_xtoy_sound(grid, iv, y)
    case ('sonde_sfc'); call da_transform_xtoy_sonde_sfc(grid, iv, y)
    case default; wrfda_xtoy_apply_grid = 1_c_int; return
    end select

    call copy_y_to_out(fam_id, var_code, y, out_y, num_obs)
    wrfda_xtoy_apply_grid = 0_c_int
  end function wrfda_xtoy_apply_grid

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

    type(domain), target :: grid
    type(iv_type), target :: iv
    type(y_type), target :: jo_grad_y
    type(x_type), target :: jo_grad_x
    integer :: n
    character(len=256) :: fambuf
    integer :: famlen
    character(len=:), allocatable :: op_str, fam_str, var_str
    integer :: fam_id
    character(len=1) :: var_code

    call init_domain_from_arrays_refs(grid, nx, ny, nz, inout_u, inout_v, inout_t, inout_q, inout_psfc)
    call zero_x_like(jo_grad_x, nx, ny, nz)
    sfc_assi_options = sfc_assi_options_1

    fambuf = '' ; famlen = 0
    do n = 1, len(fambuf)
      if (operator_family(n) == c_null_char) exit
      famlen = famlen + 1
      fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
    end do
    op_str = trim(fambuf(1:max(famlen,1)))
    call parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
    call init_iv_from_obs(fam_id, iv, nx, ny, nz, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels)
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
  integer(c_int) function wrfda_xtoy_apply_profiles(operator_family, nx, ny, nz, &
      u, v, t, q, psfc, lats2d, lons2d, levels, &
      num_obs, obs_counts, obs_lats, obs_lons, obs_levels_flat, out_y_flat) &
      bind(C, name="wrfda_xtoy_apply_profiles")
    implicit none
    character(c_char), intent(in) :: operator_family(*)
    integer(c_int), value :: nx, ny, nz
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    real(c_double), intent(in) :: lats2d(*), lons2d(*), levels(*)
    integer(c_int), value :: num_obs
    integer(c_int), intent(in) :: obs_counts(*)
    real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels_flat(*)
    real(c_double), intent(out) :: out_y_flat(*)

    type(domain), target :: grid
    type(iv_type), target :: iv
    type(y_type), target :: y
    integer :: n
    character(len=256) :: fambuf
    integer :: famlen
    character(len=:), allocatable :: op_str, fam_str, var_str
    integer :: fam_id
    character(len=1) :: var_code

    call init_domain_from_arrays(grid, nx, ny, nz, u, v, t, q, psfc)
    sfc_assi_options = sfc_assi_options_1

    fambuf = '' ; famlen = 0
    do n = 1, len(fambuf)
      if (operator_family(n) == c_null_char) exit
      famlen = famlen + 1
      fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
    end do
    op_str = trim(fambuf(1:max(famlen,1)))
    call parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
    call init_iv_profiles(fam_id, iv, nx, ny, nz, lats2d, lons2d, levels, &
                          num_obs, obs_counts, obs_lats, obs_lons, obs_levels_flat)
    call da_allocate_y(iv, y)

    select case (trim(fam_str))
    case ('airep'); call da_transform_xtoy_airep(grid, iv, y)
    case ('pilot'); call da_transform_xtoy_pilot(grid, iv, y)
    case ('sound'); call da_transform_xtoy_sound(grid, iv, y)
    case default; wrfda_xtoy_apply_profiles = 1_c_int; return
    end select

    call copy_y_to_out_profiles(fam_id, var_code, y, num_obs, obs_counts, out_y_flat)
    wrfda_xtoy_apply_profiles = 0_c_int
  end function wrfda_xtoy_apply_profiles

  integer(c_int) function wrfda_xtoy_adjoint_profiles(operator_family, nx, ny, nz, &
      delta_y_flat, lats2d, lons2d, levels, num_obs, obs_counts, obs_lats, obs_lons, obs_levels_flat, &
      inout_u, inout_v, inout_t, inout_q, inout_psfc) bind(C, name="wrfda_xtoy_adjoint_profiles")
    implicit none
    character(c_char), intent(in) :: operator_family(*)
    integer(c_int), value :: nx, ny, nz
    real(c_double), intent(in) :: delta_y_flat(*)
    real(c_double), intent(in) :: lats2d(*), lons2d(*), levels(*)
    integer(c_int), value :: num_obs
    integer(c_int), intent(in) :: obs_counts(*)
    real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels_flat(*)
    real(c_double), intent(inout) :: inout_u(*), inout_v(*), inout_t(*), inout_q(*), inout_psfc(*)

    type(domain), target :: grid
    type(iv_type), target :: iv
    type(y_type), target :: jo_grad_y
    type(x_type), target :: jo_grad_x
    integer :: n
    character(len=256) :: fambuf
    integer :: famlen
    character(len=:), allocatable :: op_str, fam_str, var_str
    integer :: fam_id
    character(len=1) :: var_code

    call init_domain_from_arrays_refs(grid, nx, ny, nz, inout_u, inout_v, inout_t, inout_q, inout_psfc)
    call zero_x_like(jo_grad_x, nx, ny, nz)
    sfc_assi_options = sfc_assi_options_1

    fambuf = '' ; famlen = 0
    do n = 1, len(fambuf)
      if (operator_family(n) == c_null_char) exit
      famlen = famlen + 1
      fambuf(famlen:famlen) = achar(iachar(operator_family(n)))
    end do
    op_str = trim(fambuf(1:max(famlen,1)))
    call parse_family_variable(op_str, fam_str, var_str, fam_id, var_code)
    call init_iv_profiles(fam_id, iv, nx, ny, nz, lats2d, lons2d, levels, &
                          num_obs, obs_counts, obs_lats, obs_lons, obs_levels_flat)
    call da_allocate_y(iv, jo_grad_y)
    call init_y_from_delta_profiles(fam_id, var_code, jo_grad_y, num_obs, obs_counts, delta_y_flat)

    select case (trim(fam_str))
    case ('airep'); call da_transform_xtoy_airep_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('pilot'); call da_transform_xtoy_pilot_adj(grid, iv, jo_grad_y, jo_grad_x)
    case ('sound'); call da_transform_xtoy_sound_adj(grid, iv, jo_grad_y, jo_grad_x)
    case default; wrfda_xtoy_adjoint_profiles = 1_c_int; return
    end select

    call copy_x_to_state(jo_grad_x, inout_u, inout_v, inout_t, inout_q, inout_psfc, nx, ny, nz)
    wrfda_xtoy_adjoint_profiles = 0_c_int
  end function wrfda_xtoy_adjoint_profiles


  subroutine init_domain_from_arrays(grid, nx, ny, nz, u, v, t, q, psfc)
    type(domain), intent(inout) :: grid
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
    ! Allocate and fill xa fields (single precision in WRFDA)
    allocate(grid%xa%u(nx,ny,nz1)); allocate(grid%xa%v(nx,ny,nz1))
    allocate(grid%xa%t(nx,ny,nz1)); allocate(grid%xa%q(nx,ny,nz1))
    allocate(grid%xa%psfc(nx,ny))
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
    allocate(grid%xb%u(nx,ny,nz1)); grid%xb%u = grid%xa%u
    allocate(grid%xb%v(nx,ny,nz1)); grid%xb%v = grid%xa%v
    allocate(grid%xb%t(nx,ny,nz1)); grid%xb%t = grid%xa%t
    allocate(grid%xb%q(nx,ny,nz1)); grid%xb%q = grid%xa%q
    allocate(grid%xb%psfc(nx,ny));  grid%xb%psfc = grid%xa%psfc
  end subroutine init_domain_from_arrays

  subroutine init_domain_from_arrays_refs(grid, nx, ny, nz, u, v, t, q, psfc)
    type(domain), intent(inout) :: grid
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(inout) :: u(*), v(*), t(*), q(*), psfc(*)
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
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
  end subroutine init_domain_from_arrays_refs

  subroutine init_iv_from_obs(family, iv, nx, ny, nz, lats2d, lons2d, levels, num_obs, obs_lats, obs_lons, obs_levels)
    integer, intent(in) :: family
    type(iv_type), intent(inout) :: iv
    integer, intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: lats2d(*), lons2d(*), levels(*)
    integer, intent(in) :: num_obs
    real(c_double), intent(in) :: obs_lats(*), obs_lons(*), obs_levels(*)
    integer :: n, i, j, k
    real(c_double) :: xfloat, yfloat, zkfloat
    ! Allocate basic fields for this family (single-level per obs)
    iv%info(family)%nlocal = num_obs
    iv%info(family)%ntotal = num_obs
    iv%info(family)%max_lev = 1
    allocate(iv%info(family)%levels(num_obs)); iv%info(family)%levels = 1
    allocate(iv%info(family)%i(1,num_obs), iv%info(family)%j(1,num_obs), iv%info(family)%k(1,num_obs))
    allocate(iv%info(family)%x(1,num_obs), iv%info(family)%y(1,num_obs), iv%info(family)%zk(1,num_obs))
    allocate(iv%info(family)%proc_domain(1,num_obs)); iv%info(family)%proc_domain = .true.
    iv%info(family)%n1 = 1; iv%info(family)%n2 = num_obs
    do n = 1, num_obs
      call find_fractional_ij(nx, ny, lats2d, lons2d, obs_lats(n), obs_lons(n), i, j, xfloat, yfloat)
      iv%info(family)%i(1,n) = i
      iv%info(family)%j(1,n) = j
      if (nz > 1) then
        call find_fractional_k(nz, levels, obs_levels(n), k, zkfloat)
      else
        k = 1; zkfloat = 1.0_c_double
      end if
      iv%info(family)%k(1,n)  = k
      iv%info(family)%x(1,n)  = real(xfloat)
      iv%info(family)%y(1,n)  = real(yfloat)
      iv%info(family)%zk(1,n) = real(zkfloat)
    end do
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
    integer :: n, i, j, k, idx, lev, nlev, max_lev
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
    allocate(iv%info(family)%levels(num_obs))
    allocate(iv%info(family)%proc_domain(max_lev,num_obs)); iv%info(family)%proc_domain = .true.
    iv%info(family)%n1 = 1; iv%info(family)%n2 = num_obs
    idx = 0
    do n = 1, num_obs
      nlev = obs_counts(n)
      iv%info(family)%levels(n) = nlev
      call find_fractional_ij(nx, ny, lats2d, lons2d, obs_lats(n), obs_lons(n), i, j, xfloat, yfloat)
      do lev = 1, nlev
        if (nz > 1) then
          call find_fractional_k(nz, levels, obs_levels_flat(idx+lev), k, zkfloat)
        else
          k = 1; zkfloat = 1.0_c_double
        end if
        iv%info(family)%i(lev,n) = i
        iv%info(family)%j(lev,n) = j
        iv%info(family)%k(lev,n) = k
        iv%info(family)%x(lev,n) = real(xfloat)
        iv%info(family)%y(lev,n) = real(yfloat)
        iv%info(family)%zk(lev,n) = real(zkfloat)
      end do
      do lev = nlev+1, max_lev
        iv%info(family)%i(lev,n) = i
        iv%info(family)%j(lev,n) = j
        iv%info(family)%k(lev,n) = 1
        iv%info(family)%x(lev,n) = real(xfloat)
        iv%info(family)%y(lev,n) = real(yfloat)
        iv%info(family)%zk(lev,n) = 1.0
      end do
      idx = idx + nlev
    end do
  end subroutine init_iv_profiles

  ! No longer used: WRFDA's da_allocate_y handles family allocation
  subroutine init_y_for_family(family, y, num_obs)
    integer, intent(in) :: family
    type(y_type), intent(inout) :: y
    integer, intent(in) :: num_obs
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
    character(len=:), allocatable, intent(out) :: fam_str, var_str
    integer, intent(out) :: fam_id
    character(len=1), intent(out) :: var_code
    integer :: pos
    fam_str = trim(op_str)
    var_str = ''
    pos = index(op_str, ':')
    if (pos > 0) then
      fam_str = trim(op_str(:pos-1))
      var_str = trim(op_str(pos+1:))
    end if
    select case (fam_str)
    case ('metar'); fam_id = metar
    case ('synop'); fam_id = synop
    case ('ships'); fam_id = ships
    case ('buoy');  fam_id = buoy
    case default;   fam_id = metar
    end select
    if (len_trim(var_str) == 0) then
      var_code = 't'
    else
      select case (var_str(1:1))
      case ('u','U'); var_code = 'u'
      case ('v','V'); var_code = 'v'
      case ('t','T'); var_code = 't'
      case ('q','Q'); var_code = 'q'
      case ('p','P'); var_code = 'p'
      case default;   var_code = 't'
      end select
    end if
  end subroutine parse_family_variable

  subroutine copy_x_to_state(jo_grad_x, u, v, t, q, psfc, nx, ny, nz)
    type(x_type), intent(in) :: jo_grad_x
    real(c_double), intent(inout) :: u(*), v(*), t(*), q(*), psfc(*)
    integer, intent(in) :: nx, ny, nz
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
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

  subroutine find_nearest_ij(nx, ny, lats2d, lons2d, olat, olon, oi, oj)
    integer, intent(in) :: nx, ny
    real(c_double), intent(in) :: lats2d(*), lons2d(*), olat, olon
    integer, intent(out) :: oi, oj
    integer :: i,j
    real(c_double) :: bestd, d, lat, lon
    bestd = 1.0d99; oi = 1; oj = 1
    do j=1,ny
      do i=1,nx
        lat = lats2d(i + (j-1)*nx)
        lon = lons2d(i + (j-1)*nx)
        d = (lat-olat)*(lat-olat) + (lon-olon)*(lon-olon)
        if (d < bestd) then
          bestd = d; oi = i; oj = j
        end if
      end do
    end do
  end subroutine find_nearest_ij

  subroutine find_nearest_k(nz, levels, olev, ok)
    integer, intent(in) :: nz
    real(c_double), intent(in) :: levels(*)
    real(c_double), intent(in) :: olev
    integer, intent(out) :: ok
    integer :: k
    real(c_double) :: bestd, d
    bestd = 1.0d99; ok = 1
    do k=1,nz
      d = abs(levels(k) - olev)
      if (d < bestd) then
        bestd = d; ok = k
      end if
    end do
  end subroutine find_nearest_k

  subroutine find_fractional_ij(nx, ny, lats2d, lons2d, olat, olon, oi, oj, xfloat, yfloat)
    integer, intent(in) :: nx, ny
    real(c_double), intent(in) :: lats2d(*), lons2d(*), olat, olon
    integer, intent(out) :: oi, oj
    real(c_double), intent(out) :: xfloat, yfloat
    integer :: i, j, inext, iprev, jnext, jprev
    real(c_double) :: bestd, d
    real(c_double) :: lat0, lon0, lat_i, lon_i, lat_j, lon_j
    real(c_double) :: dx_o, dy_o, dx_i, dy_i, dx_j, dy_j
    real(c_double) :: proj, norm2
    ! Find nearest grid point by great-circle-like proxy (equirectangular)
    bestd = 1.0d99; oi = 1; oj = 1
    do j=1,ny
      do i=1,nx
        call geo_vec(lats2d(i + (j-1)*nx), lons2d(i + (j-1)*nx), olat, olon, dx_o, dy_o)
        d = dx_o*dx_o + dy_o*dy_o
        if (d < bestd) then
          bestd = d; oi = i; oj = j
        end if
      end do
    end do
    lat0 = lats2d(oi + (oj-1)*nx); lon0 = lons2d(oi + (oj-1)*nx)
    ! Fraction along i using projection onto local i-axis
    inext = min(oi+1, nx); iprev = max(oi-1, 1)
    if (oi < nx) then
      lat_i = lats2d(inext + (oj-1)*nx); lon_i = lons2d(inext + (oj-1)*nx)
      call geo_vec(lat0, lon0, lat_i, lon_i, dx_i, dy_i)
      call geo_vec(lat0, lon0, olat,  olon,  dx_o, dy_o)
      norm2 = dx_i*dx_i + dy_i*dy_i + 1.0d-16
      proj = (dx_o*dx_i + dy_o*dy_i) / norm2
      proj = max(0.0_c_double, min(1.0_c_double, proj))
      xfloat = real(oi, c_double) + proj
    else
      lat_i = lats2d(iprev + (oj-1)*nx); lon_i = lons2d(iprev + (oj-1)*nx)
      call geo_vec(lat_i, lon_i, lat0, lon0, dx_i, dy_i)
      call geo_vec(lat_i, lon_i, olat,  olon,  dx_o, dy_o)
      norm2 = dx_i*dx_i + dy_i*dy_i + 1.0d-16
      proj = (dx_o*dx_i + dy_o*dy_i) / norm2
      proj = max(0.0_c_double, min(1.0_c_double, proj))
      xfloat = real(iprev, c_double) + proj
    end if
    ! Fraction along j using projection onto local j-axis
    jnext = min(oj+1, ny); jprev = max(oj-1, 1)
    if (oj < ny) then
      lat_j = lats2d(oi + (jnext-1)*nx); lon_j = lons2d(oi + (jnext-1)*nx)
      call geo_vec(lat0, lon0, lat_j, lon_j, dx_j, dy_j)
      call geo_vec(lat0, lon0, olat,  olon,  dx_o, dy_o)
      norm2 = dx_j*dx_j + dy_j*dy_j + 1.0d-16
      proj = (dx_o*dx_j + dy_o*dy_j) / norm2
      proj = max(0.0_c_double, min(1.0_c_double, proj))
      yfloat = real(oj, c_double) + proj
    else
      lat_j = lats2d(oi + (jprev-1)*nx); lon_j = lons2d(oi + (jprev-1)*nx)
      call geo_vec(lat_j, lon_j, lat0, lon0, dx_j, dy_j)
      call geo_vec(lat_j, lon_j, olat,  olon,  dx_o, dy_o)
      norm2 = dx_j*dx_j + dy_j*dy_j + 1.0d-16
      proj = (dx_o*dx_j + dy_o*dy_j) / norm2
      proj = max(0.0_c_double, min(1.0_c_double, proj))
      yfloat = real(jprev, c_double) + proj
    end if
  end subroutine find_fractional_ij

  subroutine find_fractional_k(nz, levels, olev, ok, zkfloat)
    integer, intent(in) :: nz
    real(c_double), intent(in) :: levels(*)
    real(c_double), intent(in) :: olev
    integer, intent(out) :: ok
    real(c_double), intent(out) :: zkfloat
    integer :: k0, k1
    real(c_double) :: d0, d1
    if (nz <= 1) then
      ok = 1; zkfloat = 1.0_c_double; return
    end if
    call find_nearest_k(nz, levels, olev, ok)
    k0 = ok
    k1 = min(k0+1, nz)
    d0 = abs(olev - levels(k0))
    d1 = abs(levels(k1) - olev) + 1.0d-12
    zkfloat = real(k0, c_double) + max(0.0_c_double, min(1.0_c_double, d0/(d0+d1)))
    if (k1 == k0) zkfloat = real(k0, c_double)
  end subroutine find_fractional_k

  subroutine geo_vec(lat0, lon0, lat1, lon1, dx, dy)
    ! Convert two lat/lon positions to local tangent-plane meters (east/north)
    real(c_double), intent(in) :: lat0, lon0, lat1, lon1
    real(c_double), intent(out) :: dx, dy
    real(c_double), parameter :: deg2rad = 0.017453292519943295d0
    real(c_double), parameter :: R = 6371000.0d0
    real(c_double) :: phi0, phi1, lam0, lam1, dphi, dlam, m
    phi0 = lat0 * deg2rad; phi1 = lat1 * deg2rad
    lam0 = lon0 * deg2rad; lam1 = lon1 * deg2rad
    dphi = phi1 - phi0
    dlam = lam1 - lam0
    m = cos(0.5d0*(phi0+phi1))
    dx = R * dlam * m
    dy = R * dphi
  end subroutine geo_vec

end module metada_wrfda_dispatch


