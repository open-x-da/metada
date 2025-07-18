!
! geometry.f90
! Object-oriented representation for the Lorenz 63 geometric properties
!

module geometry
  use state, only: state_type
  implicit none
  private
  
  ! Public type definitions and interfaces
  public :: geometry_type
  
  !> Geometry class for the Lorenz 63 system
  type :: geometry_type
    private
    real :: x_min, x_max  !< X-axis boundaries
    real :: y_min, y_max  !< Y-axis boundaries
    real :: z_min, z_max  !< Z-axis boundaries
  contains
    ! Type-bound procedures (methods)
    procedure :: init => geometry_init
    procedure :: contains_point => geometry_contains_point
    procedure :: get_x_range => geometry_get_x_range
    procedure :: get_y_range => geometry_get_y_range
    procedure :: get_z_range => geometry_get_z_range
  end type geometry_type
  
  ! Module procedures not bound to type
  public :: compute_distance
  
contains

  !> Initialize geometry bounds for the Lorenz attractor
  subroutine geometry_init(this, x_min, x_max, y_min, y_max, z_min, z_max)
    class(geometry_type), intent(inout) :: this
    real, intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
    
    this%x_min = x_min
    this%x_max = x_max
    this%y_min = y_min
    this%y_max = y_max
    this%z_min = z_min
    this%z_max = z_max
  end subroutine geometry_init
  
  !> Check if a state point is within the geometry bounds
  logical function geometry_contains_point(this, state)
    class(geometry_type), intent(in) :: this
    type(state_type), intent(in) :: state
    
    geometry_contains_point = .false.
    
    if (state%get_x() >= this%x_min .and. state%get_x() <= this%x_max .and. &
        state%get_y() >= this%y_min .and. state%get_y() <= this%y_max .and. &
        state%get_z() >= this%z_min .and. state%get_z() <= this%z_max) then
      geometry_contains_point = .true.
    end if
  end function geometry_contains_point
  
  !> Get the range of x-axis
  subroutine geometry_get_x_range(this, min_val, max_val)
    class(geometry_type), intent(in) :: this
    real, intent(out) :: min_val, max_val
    
    min_val = this%x_min
    max_val = this%x_max
  end subroutine geometry_get_x_range
  
  !> Get the range of y-axis
  subroutine geometry_get_y_range(this, min_val, max_val)
    class(geometry_type), intent(in) :: this
    real, intent(out) :: min_val, max_val
    
    min_val = this%y_min
    max_val = this%y_max
  end subroutine geometry_get_y_range
  
  !> Get the range of z-axis
  subroutine geometry_get_z_range(this, min_val, max_val)
    class(geometry_type), intent(in) :: this
    real, intent(out) :: min_val, max_val
    
    min_val = this%z_min
    max_val = this%z_max
  end subroutine geometry_get_z_range
  
  !> Compute Euclidean distance between two states in the phase space
  !> We maintain this as a module procedure for compatibility
  function compute_distance(state1, state2) result(distance)
    type(state_type), intent(in) :: state1, state2
    real :: distance
    
    ! Use the state's built-in distance method
    distance = state1%distance(state2)
  end function compute_distance
  
end module geometry 