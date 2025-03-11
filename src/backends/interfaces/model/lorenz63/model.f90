!
! model.f90
! Object-oriented implementation for the Lorenz 63 dynamical system with C bindings
!

module model
  use, intrinsic :: iso_c_binding
  use state, only: state_type
  implicit none
  private
  
  ! Public type definitions and interfaces
  public :: model_type, integrator_type, rk4_integrator_type, lorenz63_model_type
  
  ! C binding public interfaces
  public :: c_lorenz63_model_create, c_lorenz63_model_destroy
  public :: c_lorenz63_model_get_params, c_lorenz63_model_set_params
  public :: c_rk4_integrator_create, c_rk4_integrator_destroy, c_rk4_integrator_step
  
  !> Base abstract model class
  type, abstract :: model_type
    private
    real :: dt = 0.01  !< Time step for integration
  contains
    ! Type-bound procedures (methods)
    procedure(compute_derivatives_interface), deferred :: compute_derivatives
    procedure :: set_time_step => model_set_time_step
    procedure :: get_time_step => model_get_time_step
  end type model_type
  
  !> Specific Lorenz63 model implementation
  type, extends(model_type) :: lorenz63_model_type
    private
    real :: sigma = 10.0  !< Prandtl number
    real :: rho = 28.0    !< Rayleigh number
    real :: beta = 8.0/3.0 !< Geometric factor
  contains
    ! Type-bound procedures (methods)
    procedure :: init => lorenz63_init
    procedure :: compute_derivatives => lorenz63_derivatives
    procedure :: get_sigma => lorenz63_get_sigma
    procedure :: get_rho => lorenz63_get_rho
    procedure :: get_beta => lorenz63_get_beta
  end type lorenz63_model_type
  
  !> Abstract integrator class
  type, abstract :: integrator_type
  contains
    procedure(step_interface), deferred :: step
    procedure(initialize_interface), deferred :: init
  end type integrator_type
  
  !> RK4 integrator implementation
  type, extends(integrator_type) :: rk4_integrator_type
    private
    class(model_type), pointer :: model => null()
  contains
    procedure :: init => rk4_init
    procedure :: step => rk4_step
    procedure :: cleanup => rk4_cleanup
    final :: rk4_finalize
  end type rk4_integrator_type
  
  ! Abstract interface definitions
  abstract interface
    subroutine compute_derivatives_interface(this, state, dxdt, dydt, dzdt)
      import :: model_type, state_type
      class(model_type), intent(in) :: this
      type(state_type), intent(in) :: state
      real, intent(out) :: dxdt, dydt, dzdt
    end subroutine compute_derivatives_interface
    
    subroutine step_interface(this, state)
      import :: integrator_type, state_type
      class(integrator_type), intent(in) :: this
      type(state_type), intent(inout) :: state
    end subroutine step_interface
    
    subroutine initialize_interface(this, model)
      import :: integrator_type, model_type
      class(integrator_type), intent(inout) :: this
      class(model_type), target, intent(in) :: model
    end subroutine initialize_interface
  end interface
  
  ! Create aliases for backward compatibility
  public :: step
  
contains
  !> Set the time step for the model
  subroutine model_set_time_step(this, dt)
    class(model_type), intent(inout) :: this
    real, intent(in) :: dt
    
    this%dt = dt
  end subroutine model_set_time_step
  
  !> Get the time step from the model
  function model_get_time_step(this) result(dt)
    class(model_type), intent(in) :: this
    real :: dt
    
    dt = this%dt
  end function model_get_time_step
  
  !> Initialize Lorenz 63 model with given parameters
  subroutine lorenz63_init(this, sigma, rho, beta, dt)
    class(lorenz63_model_type), intent(inout) :: this
    real, intent(in) :: sigma, rho, beta, dt
    
    this%sigma = sigma
    this%rho = rho
    this%beta = beta
    call this%set_time_step(dt)
  end subroutine lorenz63_init
  
  !> Get sigma parameter
  function lorenz63_get_sigma(this) result(sigma)
    class(lorenz63_model_type), intent(in) :: this
    real :: sigma
    
    sigma = this%sigma
  end function lorenz63_get_sigma
  
  !> Get rho parameter
  function lorenz63_get_rho(this) result(rho)
    class(lorenz63_model_type), intent(in) :: this
    real :: rho
    
    rho = this%rho
  end function lorenz63_get_rho
  
  !> Get beta parameter
  function lorenz63_get_beta(this) result(beta)
    class(lorenz63_model_type), intent(in) :: this
    real :: beta
    
    beta = this%beta
  end function lorenz63_get_beta
  
  !> Compute derivatives for the Lorenz system
  subroutine lorenz63_derivatives(this, state, dxdt, dydt, dzdt)
    class(lorenz63_model_type), intent(in) :: this
    type(state_type), intent(in) :: state
    real, intent(out) :: dxdt, dydt, dzdt
    
    real :: x, y, z
    
    x = state%get_x()
    y = state%get_y()
    z = state%get_z()
    
    dxdt = this%sigma * (y - x)
    dydt = x * (this%rho - z) - y
    dzdt = x * y - this%beta * z
  end subroutine lorenz63_derivatives
  
  !> Initialize the RK4 integrator with a model
  subroutine rk4_init(this, model)
    class(rk4_integrator_type), intent(inout) :: this
    class(model_type), target, intent(in) :: model
    
    ! Clean up any existing association first
    if (associated(this%model)) then
      nullify(this%model)
    end if
    
    this%model => model
  end subroutine rk4_init
  
  !> Clean up the RK4 integrator resources
  subroutine rk4_cleanup(this)
    class(rk4_integrator_type), intent(inout) :: this
    
    if (associated(this%model)) then
      nullify(this%model)
    end if
  end subroutine rk4_cleanup
  
  !> Finalizer for the RK4 integrator
  subroutine rk4_finalize(this)
    type(rk4_integrator_type), intent(inout) :: this
    
    call this%cleanup()
  end subroutine rk4_finalize
  
  !> Advance the system one time step using 4th-order Runge-Kutta
  subroutine rk4_step(this, state)
    class(rk4_integrator_type), intent(in) :: this
    type(state_type), intent(inout) :: state
    
    real :: k1x, k1y, k1z
    real :: k2x, k2y, k2z
    real :: k3x, k3y, k3z
    real :: k4x, k4y, k4z
    type(state_type) :: temp_state
    real :: dt
    
    ! Check if model is associated
    if (.not. associated(this%model)) then
      print *, "Error: Model not initialized in RK4 integrator"
      return
    end if
    
    dt = this%model%get_time_step()
    
    ! k1 = f(y_n)
    call this%model%compute_derivatives(state, k1x, k1y, k1z)
    
    ! k2 = f(y_n + dt/2 * k1)
    call temp_state%init(state%get_x() + dt/2.0 * k1x, &
                        state%get_y() + dt/2.0 * k1y, &
                        state%get_z() + dt/2.0 * k1z)
    call this%model%compute_derivatives(temp_state, k2x, k2y, k2z)
    
    ! k3 = f(y_n + dt/2 * k2)
    call temp_state%init(state%get_x() + dt/2.0 * k2x, &
                        state%get_y() + dt/2.0 * k2y, &
                        state%get_z() + dt/2.0 * k2z)
    call this%model%compute_derivatives(temp_state, k3x, k3y, k3z)
    
    ! k4 = f(y_n + dt * k3)
    call temp_state%init(state%get_x() + dt * k3x, &
                        state%get_y() + dt * k3y, &
                        state%get_z() + dt * k3z)
    call this%model%compute_derivatives(temp_state, k4x, k4y, k4z)
    
    ! y_{n+1} = y_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    call state%set_x(state%get_x() + dt/6.0 * (k1x + 2.0*k2x + 2.0*k3x + k4x))
    call state%set_y(state%get_y() + dt/6.0 * (k1y + 2.0*k2y + 2.0*k3y + k4y))
    call state%set_z(state%get_z() + dt/6.0 * (k1z + 2.0*k2z + 2.0*k3z + k4z))
  end subroutine rk4_step
  
  !> Legacy step procedure for backward compatibility
  subroutine step(model, state)
    class(lorenz63_model_type), intent(in) :: model
    type(state_type), intent(inout) :: state
    
    type(rk4_integrator_type) :: integrator
    
    call integrator%init(model)
    call integrator%step(state)
    call integrator%cleanup()
  end subroutine step
  
  !---------------------------------------------------------------------
  ! C-Fortran interoperability bindings
  !---------------------------------------------------------------------
  
  ! Create and return a C pointer to a new lorenz63 model
  function c_lorenz63_model_create(sigma, rho, beta, dt) bind(C, name="lorenz63_model_create") result(ptr_result)
    real(c_float), value, intent(in) :: sigma, rho, beta, dt
    type(c_ptr) :: ptr_result
    type(lorenz63_model_type), pointer :: f_ptr
    
    allocate(f_ptr)
    call f_ptr%init(sigma, rho, beta, dt)
    ptr_result = c_loc(f_ptr)
  end function c_lorenz63_model_create
  
  ! Destroy a lorenz63 model created with c_lorenz63_model_create
  subroutine c_lorenz63_model_destroy(model_ptr) bind(C, name="lorenz63_model_destroy")
    type(c_ptr), value, intent(in) :: model_ptr
    type(lorenz63_model_type), pointer :: f_ptr
    
    call c_f_pointer(model_ptr, f_ptr)
    deallocate(f_ptr)
  end subroutine c_lorenz63_model_destroy
  
  ! Get the parameters of a lorenz63 model
  subroutine c_lorenz63_model_get_params(model_ptr, sigma, rho, beta, dt) bind(C, name="lorenz63_model_get_params")
    type(c_ptr), value, intent(in) :: model_ptr
    real(c_float), intent(out) :: sigma, rho, beta, dt
    type(lorenz63_model_type), pointer :: f_ptr
    
    call c_f_pointer(model_ptr, f_ptr)
    sigma = f_ptr%get_sigma()
    rho = f_ptr%get_rho()
    beta = f_ptr%get_beta()
    dt = f_ptr%get_time_step()
  end subroutine c_lorenz63_model_get_params
  
  ! Set the parameters of a lorenz63 model
  subroutine c_lorenz63_model_set_params(model_ptr, sigma, rho, beta, dt) bind(C, name="lorenz63_model_set_params")
    type(c_ptr), value, intent(in) :: model_ptr
    real(c_float), value, intent(in) :: sigma, rho, beta, dt
    type(lorenz63_model_type), pointer :: f_ptr
    
    call c_f_pointer(model_ptr, f_ptr)
    call f_ptr%init(sigma, rho, beta, dt)
  end subroutine c_lorenz63_model_set_params
  
  ! Create and return a C pointer to a new rk4 integrator
  function c_rk4_integrator_create(model_ptr) bind(C, name="rk4_integrator_create") result(ptr_result)
    type(c_ptr), value, intent(in) :: model_ptr
    type(c_ptr) :: ptr_result
    type(rk4_integrator_type), pointer :: integrator_ptr
    type(lorenz63_model_type), pointer :: model_f_ptr
    
    call c_f_pointer(model_ptr, model_f_ptr)
    
    allocate(integrator_ptr)
    call integrator_ptr%init(model_f_ptr)
    ptr_result = c_loc(integrator_ptr)
  end function c_rk4_integrator_create
  
  ! Destroy an rk4 integrator created with c_rk4_integrator_create
  subroutine c_rk4_integrator_destroy(integrator_ptr) bind(C, name="rk4_integrator_destroy")
    type(c_ptr), value, intent(in) :: integrator_ptr
    type(rk4_integrator_type), pointer :: f_ptr
    
    call c_f_pointer(integrator_ptr, f_ptr)
    call f_ptr%cleanup()
    deallocate(f_ptr)
  end subroutine c_rk4_integrator_destroy
  
  ! Perform one step of the rk4 integrator
  subroutine c_rk4_integrator_step(integrator_ptr, state_ptr) bind(C, name="rk4_integrator_step")
    type(c_ptr), value, intent(in) :: integrator_ptr, state_ptr
    type(rk4_integrator_type), pointer :: integrator_f_ptr
    type(state_type), pointer :: state_f_ptr
    
    call c_f_pointer(integrator_ptr, integrator_f_ptr)
    call c_f_pointer(state_ptr, state_f_ptr)
    
    call integrator_f_ptr%step(state_f_ptr)
  end subroutine c_rk4_integrator_step
  
end module model 