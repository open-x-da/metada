!
! geometry.f90
! Geometry representation for the Lorenz 63 chaotic system
!

module geometry
  use state, only: state_type
  implicit none
  
  private
  public :: geometry_type, init, compute_distance

  ! Type definition for geometry
  type :: geometry_type
    ! Parameters defining the phase space geometry
    real :: x_min, x_max
    real :: y_min, y_max
    real :: z_min, z_max
  end type geometry_type
  
contains

  ! Initialize geometry bounds for the Lorenz attractor
  subroutine init(geometry, x_min, x_max, y_min, y_max, z_min, z_max)
    type(geometry_type), intent(out) :: geometry
    real, intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
    
    geometry%x_min = x_min
    geometry%x_max = x_max
    geometry%y_min = y_min
    geometry%y_max = y_max
    geometry%z_min = z_min
    geometry%z_max = z_max
  end subroutine init
  
  ! Compute Euclidean distance between two states in the phase space
  function compute_distance(state1, state2) result(distance)
    type(state_type), intent(in) :: state1, state2
    real :: distance
    
    distance = sqrt((state2%x - state1%x)**2 + &
                    (state2%y - state1%y)**2 + &
                    (state2%z - state1%z)**2)
  end function compute_distance
  
end module geometry 