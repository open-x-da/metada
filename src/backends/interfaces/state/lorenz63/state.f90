!
! state.f90
! State representation for the Lorenz 63 chaotic system
!

module state
  implicit none
  
  private
  public :: state_type, init, print

  ! Type definition for state
  type :: state_type
    real :: x     ! First component of Lorenz system
    real :: y     ! Second component of Lorenz system
    real :: z     ! Third component of Lorenz system
  end type state_type
  
contains

  ! Initialize a state with given values
  subroutine init(state, x_init, y_init, z_init)
    type(state_type), intent(out) :: state
    real, intent(in) :: x_init, y_init, z_init
    
    state%x = x_init
    state%y = y_init
    state%z = z_init
  end subroutine init
  
  ! Print the current state to standard output
  subroutine print(state)
    type(state_type), intent(in) :: state
    
    print *, "Lorenz State:"
    print *, "  x = ", state%x
    print *, "  y = ", state%y
    print *, "  z = ", state%z
  end subroutine print
  
end module state 