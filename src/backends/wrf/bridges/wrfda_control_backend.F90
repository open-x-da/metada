!> @file wrfda_control_backend.F90
!> @brief Fortran/C bridges for initializing WRFDA control-variable backend
!> @details Provides minimal wrappers to set up and tear down the background
!>          error (be) structure needed for CV5 control variables without
!>          modifying native WRFDA sources.

module wrfda_control_backend_bridge
  use iso_c_binding
  use module_domain,        only : domain
  use module_configure,     only : grid_config_rec_type
  use da_define_structures, only : be_type, xbx_type,                       &
                                   da_deallocate_background_errors
  use da_setup_structures,  only : da_setup_background_errors, da_setup_cv
  use da_transfer_model,    only : da_transfer_wrftoxb
  use da_vtox_transforms,   only : da_transform_vtox, da_transform_vtox_inv, &
                                   da_transform_vtox_adj
  use da_control,           only : cv_size
  implicit none

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

    ! Reset control-variable bookkeeping before setup
    be % cv % size_jb = 0
    be % cv % size_je = 0
    be % cv % size_jp = 0
    be % cv % size_js = 0
    be % cv % size_jl = 0
    be % cv % size_jt = 0

    call da_setup_background_errors(grid, be)
    call da_setup_cv(be)

    ! Total control-vector size following da_solve logic
    be % cv % size = be % cv % size_jb + be % cv % size_je + be % cv % size_jp + &
                     be % cv % size_js + be % cv % size_jl + be % cv % size_jt

    cv_size = be % cv % size
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
  end subroutine wrfda_control_backend_finalize

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
    type(xbx_type) :: xbx
    type(grid_config_rec_type) :: config_flags
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
    do i = 1, cv_size_in
      control_single(i) = real(control(i), kind(control_single(1)))
    end do

    call da_transfer_wrftoxb(xbx, grid, config_flags)

    call da_transform_vtox(grid, cv_size_in, xbx, be, grid%ep,                &
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
    type(xbx_type) :: xbx
    type(grid_config_rec_type) :: config_flags
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
    call da_transfer_wrftoxb(xbx, grid, config_flags)
    call da_transform_vtox_inv(grid, cv_size_out, xbx, be, grid%ep,           &
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
    type(xbx_type) :: xbx
    type(grid_config_rec_type) :: config_flags
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

    call da_transfer_wrftoxb(xbx, grid, config_flags)

    call da_transform_vtox_adj(grid, cv_size_grad, xbx, be, grid%ep,          &
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
    type(xbx_type) :: xbx
    type(grid_config_rec_type) :: config_flags
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
    do i = 1, cv_size_grad
      control_single(i) = real(control_gradient(i), kind(control_single(1)))
    end do

    call da_transfer_wrftoxb(xbx, grid, config_flags)

    call da_transform_vtox(grid, cv_size_grad, xbx, be, grid%ep,              &
                           control_single, grid%vv, grid%vp)

    deallocate(control_single)
  end subroutine wrfda_control_backend_control_gradient_to_state

end module wrfda_control_backend_bridge


