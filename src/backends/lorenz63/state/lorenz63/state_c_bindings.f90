!
! state_c_bindings.f90
! C bindings for Lorenz 63 state module
!

module state_c_bindings
  use, intrinsic :: iso_c_binding
  use state, only: state_type
  implicit none
  private
  
  ! Public interfaces for C bindings
  public :: c_state_create, c_state_destroy, c_state_get_x, c_state_get_y, c_state_get_z
  public :: c_state_set_x, c_state_set_y, c_state_set_z, c_state_distance
  
contains

  !---------------------------------------------------------------------
  ! C-Fortran interoperability bindings
  !---------------------------------------------------------------------
  
  ! Create and return a C pointer to a new state
  function c_state_create(x, y, z) bind(C, name="state_create") result(ptr_result)
    real(c_float), value, intent(in) :: x, y, z
    type(c_ptr) :: ptr_result
    type(state_type), pointer :: f_ptr
    real(c_double) :: x_d, y_d, z_d
    
    allocate(f_ptr)
    x_d = real(x, c_double)
    y_d = real(y, c_double)
    z_d = real(z, c_double)
    call f_ptr%init(x_d, y_d, z_d)
    ptr_result = c_loc(f_ptr)
  end function c_state_create
  
  ! Destroy a state created with c_state_create
  subroutine c_state_destroy(state_ptr) bind(C, name="state_destroy")
    type(c_ptr), value, intent(in) :: state_ptr
    type(state_type), pointer :: f_ptr
    
    call c_f_pointer(state_ptr, f_ptr)
    deallocate(f_ptr)
  end subroutine c_state_destroy
  
  ! Get the x component of a state
  function c_state_get_x(state_ptr) bind(C, name="state_get_x") result(x)
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float) :: x
    type(state_type), pointer :: f_ptr
    real(c_double) :: x_d
    
    call c_f_pointer(state_ptr, f_ptr)
    x_d = f_ptr%get_x()
    x = real(x_d, c_float)
  end function c_state_get_x
  
  ! Get the y component of a state
  function c_state_get_y(state_ptr) bind(C, name="state_get_y") result(y)
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float) :: y
    type(state_type), pointer :: f_ptr
    real(c_double) :: y_d
    
    call c_f_pointer(state_ptr, f_ptr)
    y_d = f_ptr%get_y()
    y = real(y_d, c_float)
  end function c_state_get_y
  
  ! Get the z component of a state
  function c_state_get_z(state_ptr) bind(C, name="state_get_z") result(z)
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float) :: z
    type(state_type), pointer :: f_ptr
    real(c_double) :: z_d
    
    call c_f_pointer(state_ptr, f_ptr)
    z_d = f_ptr%get_z()
    z = real(z_d, c_float)
  end function c_state_get_z
  
  ! Set the x component of a state
  subroutine c_state_set_x(state_ptr, x) bind(C, name="state_set_x")
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float), value, intent(in) :: x
    type(state_type), pointer :: f_ptr
    real(c_double) :: x_d
    
    call c_f_pointer(state_ptr, f_ptr)
    x_d = real(x, c_double)
    call f_ptr%set_x(x_d)
  end subroutine c_state_set_x
  
  ! Set the y component of a state
  subroutine c_state_set_y(state_ptr, y) bind(C, name="state_set_y")
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float), value, intent(in) :: y
    type(state_type), pointer :: f_ptr
    real(c_double) :: y_d
    
    call c_f_pointer(state_ptr, f_ptr)
    y_d = real(y, c_double)
    call f_ptr%set_y(y_d)
  end subroutine c_state_set_y
  
  ! Set the z component of a state
  subroutine c_state_set_z(state_ptr, z) bind(C, name="state_set_z")
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float), value, intent(in) :: z
    type(state_type), pointer :: f_ptr
    real(c_double) :: z_d
    
    call c_f_pointer(state_ptr, f_ptr)
    z_d = real(z, c_double)
    call f_ptr%set_z(z_d)
  end subroutine c_state_set_z
  
  ! Compute distance between two states
  function c_state_distance(state_ptr1, state_ptr2) bind(C, name="state_distance") result(dist)
    type(c_ptr), value, intent(in) :: state_ptr1, state_ptr2
    real(c_float) :: dist
    type(state_type), pointer :: f_ptr1, f_ptr2
    real(c_double) :: dist_d
    
    call c_f_pointer(state_ptr1, f_ptr1)
    call c_f_pointer(state_ptr2, f_ptr2)
    dist_d = f_ptr1%distance(f_ptr2)
    dist = real(dist_d, c_float)
  end function c_state_distance
  
end module state_c_bindings 