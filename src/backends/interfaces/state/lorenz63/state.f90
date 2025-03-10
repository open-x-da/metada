!
! state.f90
! Object-oriented representation for the Lorenz 63 chaotic system state
!

module state
  implicit none
  private
  
  ! Public type definitions and interfaces
  public :: state_type

  !> State class for the Lorenz 63 system
  type :: state_type
    private
    real :: x     !< First component of Lorenz system
    real :: y     !< Second component of Lorenz system 
    real :: z     !< Third component of Lorenz system
  contains
    ! Type-bound procedures (methods)
    procedure :: init => state_init
    procedure :: print => state_print
    procedure :: get_x => state_get_x
    procedure :: get_y => state_get_y
    procedure :: get_z => state_get_z
    procedure :: set_x => state_set_x
    procedure :: set_y => state_set_y
    procedure :: set_z => state_set_z
    procedure :: distance => state_distance
    generic :: operator(.distance.) => distance
    procedure :: copy => state_copy
  end type state_type
  
  ! Define proper interface for assignment
  interface assignment(=)
    module procedure state_assign
  end interface
  
contains

  !> Initialize a state with given values
  subroutine state_init(this, x_init, y_init, z_init)
    class(state_type), intent(inout) :: this
    real, intent(in) :: x_init, y_init, z_init
    
    this%x = x_init
    this%y = y_init
    this%z = z_init
  end subroutine state_init
  
  !> Print the current state to standard output
  subroutine state_print(this)
    class(state_type), intent(in) :: this
    
    print *, "Lorenz State:"
    print *, "  x = ", this%x
    print *, "  y = ", this%y
    print *, "  z = ", this%z
  end subroutine state_print
  
  !> Get methods for state components
  function state_get_x(this) result(x)
    class(state_type), intent(in) :: this
    real :: x
    x = this%x
  end function state_get_x
  
  function state_get_y(this) result(y)
    class(state_type), intent(in) :: this
    real :: y
    y = this%y
  end function state_get_y
  
  function state_get_z(this) result(z)
    class(state_type), intent(in) :: this
    real :: z
    z = this%z
  end function state_get_z
  
  !> Set methods for state components
  subroutine state_set_x(this, x)
    class(state_type), intent(inout) :: this
    real, intent(in) :: x
    this%x = x
  end subroutine state_set_x
  
  subroutine state_set_y(this, y)
    class(state_type), intent(inout) :: this
    real, intent(in) :: y
    this%y = y
  end subroutine state_set_y
  
  subroutine state_set_z(this, z)
    class(state_type), intent(inout) :: this
    real, intent(in) :: z
    this%z = z
  end subroutine state_set_z
  
  !> Compute Euclidean distance between two states
  real function state_distance(this, other)
    class(state_type), intent(in) :: this
    class(state_type), intent(in) :: other
    
    state_distance = sqrt((other%x - this%x)**2 + &
                         (other%y - this%y)**2 + &
                         (other%z - this%z)**2)
  end function state_distance
  
  !> Copy method as an alternative to assignment
  subroutine state_copy(this, other)
    class(state_type), intent(inout) :: this
    class(state_type), intent(in) :: other
    
    this%x = other%x
    this%y = other%y
    this%z = other%z
  end subroutine state_copy
  
  !> Assignment operator (module procedure, not type-bound)
  subroutine state_assign(lhs, rhs)
    type(state_type), intent(out) :: lhs
    type(state_type), intent(in) :: rhs
    
    lhs%x = rhs%x
    lhs%y = rhs%y
    lhs%z = rhs%z
  end subroutine state_assign
  
end module state 