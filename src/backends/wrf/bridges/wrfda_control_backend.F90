!> @file wrfda_control_backend.F90
!> @brief Fortran/C bridges for initializing WRFDA control-variable backend
!> @details Provides minimal wrappers to set up and tear down the background
!>          error (be) structure needed for CV5 control variables without
!>          modifying native WRFDA sources.

module wrfda_control_backend_bridge
  use iso_c_binding
  use module_domain,        only : domain
  use module_dm,            only : wrf_dm_sum_real, wrf_dm_sum_integer
  use da_define_structures, only : be_type, xbx_type,                       &
                                   da_deallocate_background_errors,        &
                                   da_initialize_cv, da_zero_vp_type
  use da_setup_structures,  only : da_setup_background_errors
  use da_vtox_transforms,   only : da_transform_vtox, da_transform_vtox_inv, &
                                   da_transform_vtox_adj
  use da_control,           only : cv_size, trace_use, test_dm_exact, &
                                   rootproc, cv_size_domain
  use da_tools,             only : da_trace_entry, da_trace_exit
  use da_par_util,          only : da_cv_to_global
  use da_wrf_interfaces,    only : wrf_dm_bcast_real
  implicit none

  real, allocatable, save :: wrfda_control_scratch(:)
  type(xbx_type), pointer, save :: backend_xbx => null()
  logical, save :: backend_xbx_initialized = .false.
  type(be_type), pointer, save :: backend_be => null()

  interface
    function wrfda_get_persistent_xbx() result(xbx_ptr) bind(C, name="wrfda_get_persistent_xbx")
      use iso_c_binding
      type(c_ptr) :: xbx_ptr
    end function wrfda_get_persistent_xbx
  end interface

contains

  !> @brief Initialize background-error structure for control-variable backend
  subroutine wrfda_control_backend_setup(grid_ptr, be_ptr, cv_size_out, error_code) &
      bind(C, name="wrfda_control_backend_setup")
    type(c_ptr), value       :: grid_ptr
    type(c_ptr), intent(out) :: be_ptr
    integer(c_int), intent(out) :: cv_size_out
    integer(c_int), intent(out) :: error_code

    type(domain), pointer :: grid
    type(be_type), pointer :: be

    error_code   = 0
    be_ptr       = c_null_ptr
    cv_size_out  = 0

    if (.not. c_associated(grid_ptr)) then
      error_code = 1
      return
    end if

    call c_f_pointer(grid_ptr, grid)

    allocate(be)
    backend_be => be

    ! Reset control-variable bookkeeping before setup
    be % cv % size_jb = 0
    be % cv % size_je = 0
    be % cv % size_jp = 0
    be % cv % size_js = 0
    be % cv % size_jl = 0
    be % cv % size_jt = 0

    call da_setup_background_errors(grid, be)
    ! Note: da_setup_background_errors internally calls da_setup_cv(be) at line 117
    ! So we don't need to call it separately here

    ! Total control-vector size following da_solve logic
    ! Compute cv_size early so we can use it for initialization
    be % cv % size = be % cv % size_jb + be % cv % size_je + be % cv % size_jp + &
                     be % cv % size_js + be % cv % size_jl + be % cv % size_jt

    cv_size = be % cv % size

    ! Match WRFDA's initialization sequence from da_solve.inc:
    ! 1. Initialize control variable (cvt) before zeroing work arrays
    ! 2. Zero work arrays (vp, vv)
    ! 3. Initialize another control variable (xhat) after zeroing
    ! This matches the exact sequence in WRFDA trace (lines 2541-2548)
    if (cv_size > 0) then
      if (allocated(wrfda_control_scratch)) then
        if (size(wrfda_control_scratch) /= cv_size) then
          deallocate(wrfda_control_scratch)
        end if
      end if
      if (.not. allocated(wrfda_control_scratch)) then
        allocate(wrfda_control_scratch(cv_size))
      end if
      ! First da_initialize_cv call (for cvt-equivalent) before da_zero_vp_type
      ! This matches WRFDA's da_initialize_cv(cv_size, cvt) at line 655
      call da_initialize_cv(cv_size, wrfda_control_scratch)
    end if

    ! Zero domain control-variable work arrays
    ! Note: vp and vv are domain-level work arrays used by control transforms
    ! This matches WRFDA's da_zero_vp_type calls at lines 656-657
    call da_zero_vp_type(grid%vp)
    call da_zero_vp_type(grid%vv)

    ! Second da_initialize_cv call (for xhat-equivalent) after da_zero_vp_type
    ! In WRFDA, this initializes xhat which is used in the minimization loop
    ! This matches WRFDA's da_initialize_cv(cv_size, xhat) at line 805
    if (cv_size > 0 .and. allocated(wrfda_control_scratch)) then
      call da_initialize_cv(cv_size, wrfda_control_scratch)
    end if

    call ensure_backend_xbx_ready()
    cv_size_out = cv_size
    be_ptr      = c_loc(be)
  end subroutine wrfda_control_backend_setup

  !> @brief Release background-error structure allocated by setup
  subroutine wrfda_control_backend_finalize(be_ptr, error_code) &
      bind(C, name="wrfda_control_backend_finalize")
    type(c_ptr), value :: be_ptr
    integer(c_int), intent(out) :: error_code

    type(be_type), pointer :: be

    error_code = 0

    if (.not. c_associated(be_ptr)) return

    call c_f_pointer(be_ptr, be)

    if (associated(be)) then
      call da_deallocate_background_errors(be)
      deallocate(be)
    end if

    if (allocated(wrfda_control_scratch)) then
      deallocate(wrfda_control_scratch)
    end if

    backend_xbx_initialized = .false.
    if (associated(backend_xbx)) then
      nullify(backend_xbx)
    end if
    if (associated(backend_be)) then
      nullify(backend_be)
    end if
  end subroutine wrfda_control_backend_finalize

  !> @brief Return BE pointer saved during setup (for Fortran-side callers)
  function wrfda_control_backend_get_be_ptr() result(be_cptr) bind(C, name="wrfda_control_backend_get_be_ptr")
    use iso_c_binding, only : c_ptr, c_loc, c_null_ptr
    type(c_ptr) :: be_cptr
    if (associated(backend_be)) then
      be_cptr = c_loc(backend_be)
    else
      be_cptr = c_null_ptr
    end if
  end function wrfda_control_backend_get_be_ptr

  !> @brief Return current control vector size (CV5) from backend BE
  integer(c_int) function wrfda_control_backend_get_cv_size() bind(C, name="wrfda_control_backend_get_cv_size")
    use iso_c_binding, only : c_int
    use da_control, only : cv_size
    if (associated(backend_be)) then
      wrfda_control_backend_get_cv_size = backend_be % cv % size
    else
      wrfda_control_backend_get_cv_size = cv_size
    end if
  end function wrfda_control_backend_get_cv_size

  subroutine wrfda_control_backend_control_to_state(grid_ptr, be_ptr, control, &
                                                   cv_size_in, error_code)      &
      bind(C, name="wrfda_control_backend_control_to_state")
    use iso_c_binding, only : c_ptr, c_f_pointer, c_associated, c_int, c_double
    use module_domain, only : domain
    use da_define_structures, only : be_type
    implicit none

    type(c_ptr), value :: grid_ptr
    type(c_ptr), value :: be_ptr
    real(c_double), intent(in) :: control(*)
    integer(c_int), value :: cv_size_in
    integer(c_int), intent(out) :: error_code

    type(domain), pointer :: grid
    type(be_type), pointer :: be
    real, allocatable :: control_single(:)
    integer :: i

    error_code = 0

    if (.not. c_associated(grid_ptr)) then
      error_code = 1
      return
    end if

    if (.not. c_associated(be_ptr)) then
      error_code = 2
      return
    end if

    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(be_ptr, be)

    if (cv_size_in <= 0) then
      error_code = 3
      return
    end if

    allocate(control_single(cv_size_in))
    control_single = 0.0
    do i = 1, cv_size_in
      control_single(i) = real(control(i), kind(control_single(1)))
    end do

    call da_transform_vtox(grid, cv_size_in, backend_xbx, be, grid%ep,        &
                           control_single, grid%vv, grid%vp)

    deallocate(control_single)
  end subroutine wrfda_control_backend_control_to_state

  subroutine wrfda_control_backend_state_to_control(grid_ptr, be_ptr, control, &
                                                   cv_size_out, error_code)     &
      bind(C, name="wrfda_control_backend_state_to_control")
    use iso_c_binding, only : c_ptr, c_f_pointer, c_associated, c_int, c_double
    use module_domain, only : domain
    use da_define_structures, only : be_type
    implicit none

    type(c_ptr), value :: grid_ptr
    type(c_ptr), value :: be_ptr
    real(c_double), intent(out) :: control(*)
    integer(c_int), value :: cv_size_out
    integer(c_int), intent(out) :: error_code

    type(domain), pointer :: grid
    type(be_type), pointer :: be
    real, allocatable :: control_single(:)
    integer :: i

    error_code = 0

    if (.not. c_associated(grid_ptr)) then
      error_code = 1
      return
    end if

    if (.not. c_associated(be_ptr)) then
      error_code = 2
      return
    end if

    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(be_ptr, be)

    if (cv_size_out <= 0) then
      error_code = 3
      return
    end if

    allocate(control_single(cv_size_out))
    control_single = 0.0
    call da_transform_vtox_inv(grid, cv_size_out, backend_xbx, be, grid%ep,   &
                               control_single, grid%vv, grid%vp)

    do i = 1, cv_size_out
      control(i) = control_single(i)
    end do

    deallocate(control_single)
  end subroutine wrfda_control_backend_state_to_control

  subroutine wrfda_control_backend_state_gradient_to_control(grid_ptr, be_ptr, &
                                                            control_gradient, &
                                                            cv_size_grad,     &
                                                            error_code)       &
      bind(C, name="wrfda_control_backend_state_gradient_to_control")
    use iso_c_binding, only : c_ptr, c_f_pointer, c_associated, c_int, c_double
    use module_domain, only : domain
    use da_define_structures, only : be_type
    implicit none

    type(c_ptr), value :: grid_ptr
    type(c_ptr), value :: be_ptr
    real(c_double), intent(out) :: control_gradient(*)
    integer(c_int), value :: cv_size_grad
    integer(c_int), intent(out) :: error_code

    type(domain), pointer :: grid
    type(be_type), pointer :: be
    real, allocatable :: control_single(:)
    integer :: i

    error_code = 0

    if (.not. c_associated(grid_ptr)) then
      error_code = 1
      return
    end if

    if (.not. c_associated(be_ptr)) then
      error_code = 2
      return
    end if

    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(be_ptr, be)

    if (cv_size_grad <= 0) then
      error_code = 3
      return
    end if

    allocate(control_single(cv_size_grad))
    control_single = 0.0
    call da_transform_vtox_adj(grid, cv_size_grad, backend_xbx, be, grid%ep,  &
                               grid%vp, grid%vv, control_single)

    do i = 1, cv_size_grad
      control_gradient(i) = control_single(i)
    end do

    deallocate(control_single)
  end subroutine wrfda_control_backend_state_gradient_to_control

  subroutine wrfda_control_backend_control_gradient_to_state(grid_ptr, be_ptr, &
                                                            control_gradient, &
                                                            cv_size_grad,     &
                                                            error_code)       &
      bind(C, name="wrfda_control_backend_control_gradient_to_state")
    use iso_c_binding, only : c_ptr, c_f_pointer, c_associated, c_int, c_double
    use module_domain, only : domain
    use da_define_structures, only : be_type
    implicit none

    type(c_ptr), value :: grid_ptr
    type(c_ptr), value :: be_ptr
    real(c_double), intent(in) :: control_gradient(*)
    integer(c_int), value :: cv_size_grad
    integer(c_int), intent(out) :: error_code

    type(domain), pointer :: grid
    type(be_type), pointer :: be
    real, allocatable :: control_single(:)
    integer :: i

    error_code = 0

    if (.not. c_associated(grid_ptr)) then
      error_code = 1
      return
    end if

    if (.not. c_associated(be_ptr)) then
      error_code = 2
      return
    end if

    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(be_ptr, be)

    if (cv_size_grad <= 0) then
      error_code = 3
      return
    end if

    allocate(control_single(cv_size_grad))
    control_single = 0.0
    do i = 1, cv_size_grad
      control_single(i) = real(control_gradient(i), kind(control_single(1)))
    end do

    call da_transform_vtox(grid, cv_size_grad, backend_xbx, be, grid%ep,      &
                           control_single, grid%vv, grid%vp)

    deallocate(control_single)
  end subroutine wrfda_control_backend_control_gradient_to_state

  subroutine wrfda_control_backend_control_dot(grid_ptr, be_ptr, control_a,   &
                                               control_b, cv_size_in,         &
                                               dot_value, error_code)         &
      bind(C, name="wrfda_control_backend_control_dot")
    use iso_c_binding, only : c_ptr, c_f_pointer, c_associated, c_int,        &
                              c_double
    use module_domain, only : domain
    use module_dm, only : wrf_dm_sum_real
    use da_define_structures, only : be_type
    implicit none

    type(c_ptr), value :: grid_ptr
    type(c_ptr), value :: be_ptr
    real(c_double), intent(in) :: control_a(*)
    real(c_double), intent(in) :: control_b(*)
    integer(c_int), value :: cv_size_in
    real(c_double), intent(out) :: dot_value
    integer(c_int), intent(out) :: error_code

    type(domain), pointer :: grid
    type(be_type), pointer :: be
    real, allocatable :: control_single_a(:)
    real, allocatable :: control_single_b(:)
    integer :: i
    real :: dot_single

    error_code = 0
    dot_value = 0.0_c_double

    if (.not. c_associated(grid_ptr)) then
      error_code = 1
      return
    end if

    if (.not. c_associated(be_ptr)) then
      error_code = 2
      return
    end if

    if (cv_size_in <= 0) then
      error_code = 3
      return
    end if

    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(be_ptr, be)

    allocate(control_single_a(cv_size_in))
    allocate(control_single_b(cv_size_in))
    control_single_a = 0.0
    control_single_b = 0.0

    do i = 1, cv_size_in
      control_single_a(i) = real(control_a(i), kind(control_single_a(1)))
      control_single_b(i) = real(control_b(i), kind(control_single_b(1)))
    end do

    ! Use WRFDA's da_dot_cv function for proper control-space dot product
    ! This function is included from da_dot_cv.inc as a contained procedure
    dot_single = da_dot_cv(cv_size_in, control_single_a, control_single_b, &
                           grid, be%cv_mz, be%ncv_mz)
    dot_value = real(dot_single, kind=c_double)

    deallocate(control_single_a)
    deallocate(control_single_b)
  end subroutine wrfda_control_backend_control_dot

  subroutine ensure_backend_xbx_ready()
    type(c_ptr) :: xbx_ptr

    if (backend_xbx_initialized) then
      return
    end if

    xbx_ptr = wrfda_get_persistent_xbx()
    if (.not. c_associated(xbx_ptr)) then
      return
    end if

    call c_f_pointer(xbx_ptr, backend_xbx)
    backend_xbx_initialized = associated(backend_xbx)
  end subroutine ensure_backend_xbx_ready

  ! Include da_dot.inc first (needed by da_dot_cv)
  ! This defines da_dot function
  ! Path: from src/backends/wrf/bridges/ go up 5 levels to workspace root, then to WRF
#include "../../../../../WRF/var/da/da_minimisation/da_dot.inc"
  
  ! Include da_dot_cv.inc (the actual implementation)
  ! This defines da_dot_cv function
#include "../../../../../WRF/var/da/da_minimisation/da_dot_cv.inc"

end module wrfda_control_backend_bridge


