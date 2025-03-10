!
! model.f90
! Model implementation for the Lorenz 63 chaotic system
!

module model
  use state, only: state_type
  implicit none
  
  private
  public :: model_type, init, step

  ! Type definition for model parameters
  type :: model_type
    real :: sigma  ! Prandtl number
    real :: rho    ! Rayleigh number
    real :: beta   ! Geometric factor
    real :: dt     ! Time step for integration
  end type model_type
  
contains

  ! Initialize model parameters
  subroutine init(model, sigma, rho, beta, dt)
    type(model_type), intent(out) :: model
    real, intent(in) :: sigma, rho, beta, dt
    
    model%sigma = sigma
    model%rho = rho
    model%beta = beta
    model%dt = dt
  end subroutine init
  
  ! Compute derivatives for the Lorenz system
  subroutine derivatives(model, state, dxdt, dydt, dzdt)
    type(model_type), intent(in) :: model
    type(state_type), intent(in) :: state
    real, intent(out) :: dxdt, dydt, dzdt
    
    dxdt = model%sigma * (state%y - state%x)
    dydt = state%x * (model%rho - state%z) - state%y
    dzdt = state%x * state%y - model%beta * state%z
  end subroutine derivatives
  
  ! Advance the system one time step using 4th-order Runge-Kutta
  subroutine step(model, state)
    type(model_type), intent(in) :: model
    type(state_type), intent(inout) :: state
    
    real :: k1x, k1y, k1z
    real :: k2x, k2y, k2z
    real :: k3x, k3y, k3z
    real :: k4x, k4y, k4z
    type(state_type) :: temp_state
    
    ! k1 = f(y_n)
    call derivatives(model, state, k1x, k1y, k1z)
    
    ! k2 = f(y_n + dt/2 * k1)
    temp_state%x = state%x + model%dt/2.0 * k1x
    temp_state%y = state%y + model%dt/2.0 * k1y
    temp_state%z = state%z + model%dt/2.0 * k1z
    call derivatives(model, temp_state, k2x, k2y, k2z)
    
    ! k3 = f(y_n + dt/2 * k2)
    temp_state%x = state%x + model%dt/2.0 * k2x
    temp_state%y = state%y + model%dt/2.0 * k2y
    temp_state%z = state%z + model%dt/2.0 * k2z
    call derivatives(model, temp_state, k3x, k3y, k3z)
    
    ! k4 = f(y_n + dt * k3)
    temp_state%x = state%x + model%dt * k3x
    temp_state%y = state%y + model%dt * k3y
    temp_state%z = state%z + model%dt * k3z
    call derivatives(model, temp_state, k4x, k4y, k4z)
    
    ! y_{n+1} = y_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    state%x = state%x + model%dt/6.0 * (k1x + 2.0*k2x + 2.0*k3x + k4x)
    state%y = state%y + model%dt/6.0 * (k1y + 2.0*k2y + 2.0*k3y + k4y)
    state%z = state%z + model%dt/6.0 * (k1z + 2.0*k2z + 2.0*k3z + k4z)
  end subroutine step
  
end module model 