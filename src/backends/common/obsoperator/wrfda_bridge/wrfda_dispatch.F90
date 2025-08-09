! Dispatch C-bindings that translate flat arrays into WRFDA types
! and call the included transform routines.

module metada_wrfda_dispatch
  use iso_c_binding
  use module_domain,        only: domain, x_type
  use da_define_structures, only: iv_type, y_type
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
    ! Stub implementation: mark dummies as used and set out_y to zero to avoid warnings
    if (.false.) then
      print *, operator_family(1), wrfda_root(1)
      print *, state_values(1), obs_lats(1), obs_lons(1), obs_levels(1)
      print *, nx, ny, nz, num_obs
    end if
    if (num_obs > 0_c_int) out_y(1:num_obs) = 0.0_c_double
    wrfda_xtoy_apply = 1_c_int
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
    ! Stub implementation: mark dummies as used to avoid warnings
    if (.false.) then
      print *, operator_family(1), wrfda_root(1)
      print *, delta_y(1), obs_lats(1), obs_lons(1), obs_levels(1)
      print *, num_obs, nx, ny, nz, inout_state_values(1)
    end if
    wrfda_xtoy_adjoint = 1_c_int
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

end module metada_wrfda_dispatch


