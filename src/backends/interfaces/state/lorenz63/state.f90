!
! state.f90
! Object-oriented representation for the Lorenz 63 chaotic system state with C bindings
!

module state
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  ! Public type definitions and interfaces
  public :: state_type
  public :: c_state_create, c_state_destroy, c_state_get_x, c_state_get_y, c_state_get_z
  public :: c_state_set_x, c_state_set_y, c_state_set_z, c_state_distance

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

  ! Make assignment operator public  
  public :: assignment(=)
  
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
  
  !---------------------------------------------------------------------
  ! C-Fortran interoperability bindings
  !---------------------------------------------------------------------
  
  ! Create and return a C pointer to a new state
  function c_state_create(x, y, z) bind(C, name="state_create") result(ptr_result)
    real(c_float), value, intent(in) :: x, y, z
    type(c_ptr) :: ptr_result
    type(state_type), pointer :: f_ptr
    
    allocate(f_ptr)
    call f_ptr%init(x, y, z)
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
    
    call c_f_pointer(state_ptr, f_ptr)
    x = f_ptr%get_x()
  end function c_state_get_x
  
  ! Get the y component of a state
  function c_state_get_y(state_ptr) bind(C, name="state_get_y") result(y)
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float) :: y
    type(state_type), pointer :: f_ptr
    
    call c_f_pointer(state_ptr, f_ptr)
    y = f_ptr%get_y()
  end function c_state_get_y
  
  ! Get the z component of a state
  function c_state_get_z(state_ptr) bind(C, name="state_get_z") result(z)
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float) :: z
    type(state_type), pointer :: f_ptr
    
    call c_f_pointer(state_ptr, f_ptr)
    z = f_ptr%get_z()
  end function c_state_get_z
  
  ! Set the x component of a state
  subroutine c_state_set_x(state_ptr, x) bind(C, name="state_set_x")
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float), value, intent(in) :: x
    type(state_type), pointer :: f_ptr
    
    call c_f_pointer(state_ptr, f_ptr)
    call f_ptr%set_x(x)
  end subroutine c_state_set_x
  
  ! Set the y component of a state
  subroutine c_state_set_y(state_ptr, y) bind(C, name="state_set_y")
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float), value, intent(in) :: y
    type(state_type), pointer :: f_ptr
    
    call c_f_pointer(state_ptr, f_ptr)
    call f_ptr%set_y(y)
  end subroutine c_state_set_y
  
  ! Set the z component of a state
  subroutine c_state_set_z(state_ptr, z) bind(C, name="state_set_z")
    type(c_ptr), value, intent(in) :: state_ptr
    real(c_float), value, intent(in) :: z
    type(state_type), pointer :: f_ptr
    
    call c_f_pointer(state_ptr, f_ptr)
    call f_ptr%set_z(z)
  end subroutine c_state_set_z
  
  ! Compute distance between two states
  function c_state_distance(state_ptr1, state_ptr2) bind(C, name="state_distance") result(dist)
    type(c_ptr), value, intent(in) :: state_ptr1, state_ptr2
    real(c_float) :: dist
    type(state_type), pointer :: f_ptr1, f_ptr2
    
    call c_f_pointer(state_ptr1, f_ptr1)
    call c_f_pointer(state_ptr2, f_ptr2)
    dist = f_ptr1%distance(f_ptr2)
  end function c_state_distance
  
end module state 