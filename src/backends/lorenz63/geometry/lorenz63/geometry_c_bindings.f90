!
! geometry_c_bindings.f90
! C bindings for Lorenz 63 geometry module
!

module geometry_c_bindings
  use, intrinsic :: iso_c_binding
  use geometry, only: geometry_type, compute_distance
  use state, only: state_type
  implicit none
  private
  
  ! Public interfaces for C bindings
  public :: c_geometry_create, c_geometry_destroy
  public :: c_geometry_contains_point, c_geometry_get_x_range
  public :: c_geometry_get_y_range, c_geometry_get_z_range
  
contains

  !> Create a new geometry object with specified bounds
  function c_geometry_create(x_min, x_max, y_min, y_max, z_min, z_max) bind(C, name="geometry_create") result(ptr)
    real(c_float), value, intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
    type(c_ptr) :: ptr
    
    type(geometry_type), pointer :: geom
    real :: x_min_f, x_max_f, y_min_f, y_max_f, z_min_f, z_max_f
    
    ! Allocate new geometry object
    allocate(geom)
    
    ! Convert C float inputs to Fortran real (single precision)
    x_min_f = real(x_min)
    x_max_f = real(x_max)
    y_min_f = real(y_min)
    y_max_f = real(y_max)
    z_min_f = real(z_min)
    z_max_f = real(z_max)
    
    ! Initialize with provided bounds
    call geom%init(x_min_f, x_max_f, y_min_f, y_max_f, z_min_f, z_max_f)
    
    ! Return pointer to C
    ptr = c_loc(geom)
  end function c_geometry_create
  
  !> Destroy a geometry object
  subroutine c_geometry_destroy(geom_ptr) bind(C, name="geometry_destroy")
    type(c_ptr), value, intent(in) :: geom_ptr
    
    type(geometry_type), pointer :: geom
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(geom_ptr, geom)
    
    ! Deallocate the geometry object
    deallocate(geom)
  end subroutine c_geometry_destroy
  
  !> Check if a state is within the geometry bounds
  function c_geometry_contains_point(geom_ptr, state_ptr) bind(C, name="geometry_contains_point") result(result)
    type(c_ptr), value, intent(in) :: geom_ptr, state_ptr
    integer(c_int) :: result
    
    type(geometry_type), pointer :: geom
    type(state_type), pointer :: state
    
    ! Convert C pointers to Fortran pointers
    call c_f_pointer(geom_ptr, geom)
    call c_f_pointer(state_ptr, state)
    
    ! Check if state is within bounds
    if (geom%contains_point(state)) then
      result = 1
    else
      result = 0
    end if
  end function c_geometry_contains_point
  
  !> Get the x-axis range
  subroutine c_geometry_get_x_range(geom_ptr, min_val, max_val) bind(C, name="geometry_get_x_range")
    type(c_ptr), value, intent(in) :: geom_ptr
    real(c_float), intent(out) :: min_val, max_val
    
    type(geometry_type), pointer :: geom
    real :: min_val_f, max_val_f
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(geom_ptr, geom)
    
    ! Get the range
    call geom%get_x_range(min_val_f, max_val_f)
    
    ! Convert Fortran real outputs to C float
    min_val = real(min_val_f, c_float)
    max_val = real(max_val_f, c_float)
  end subroutine c_geometry_get_x_range
  
  !> Get the y-axis range
  subroutine c_geometry_get_y_range(geom_ptr, min_val, max_val) bind(C, name="geometry_get_y_range")
    type(c_ptr), value, intent(in) :: geom_ptr
    real(c_float), intent(out) :: min_val, max_val
    
    type(geometry_type), pointer :: geom
    real :: min_val_f, max_val_f
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(geom_ptr, geom)
    
    ! Get the range
    call geom%get_y_range(min_val_f, max_val_f)
    
    ! Convert Fortran real outputs to C float
    min_val = real(min_val_f, c_float)
    max_val = real(max_val_f, c_float)
  end subroutine c_geometry_get_y_range
  
  !> Get the z-axis range
  subroutine c_geometry_get_z_range(geom_ptr, min_val, max_val) bind(C, name="geometry_get_z_range")
    type(c_ptr), value, intent(in) :: geom_ptr
    real(c_float), intent(out) :: min_val, max_val
    
    type(geometry_type), pointer :: geom
    real :: min_val_f, max_val_f
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(geom_ptr, geom)
    
    ! Get the range
    call geom%get_z_range(min_val_f, max_val_f)
    
    ! Convert Fortran real outputs to C float
    min_val = real(min_val_f, c_float)
    max_val = real(max_val_f, c_float)
  end subroutine c_geometry_get_z_range
  
end module geometry_c_bindings 