!
! main.f90
! Main program to demonstrate the Lorenz 63 chaotic system using OOP
!

program main
  use state, only: state_type
  use geometry, only: geometry_type, compute_distance
  use model, only: lorenz63_model_type, rk4_integrator_type
  implicit none
  
  ! Create instances of our classes
  type(state_type) :: current_state, initial_state
  type(geometry_type) :: geom
  type(lorenz63_model_type) :: lorenz_model
  type(rk4_integrator_type) :: integrator
  
  integer :: i, num_steps
  real :: distance
  
  ! Initialize the system
  call current_state%init(1.0, 1.0, 1.0)
  call initial_state%init(1.0, 1.0, 1.0) ! Keep a copy of initial state
  
  ! Standard parameters for Lorenz 63 (sigma=10, rho=28, beta=8/3)
  call lorenz_model%init(10.0, 28.0, 8.0/3.0, 0.01)
  
  ! Initialize the integrator with our model
  call integrator%init(lorenz_model)
  
  ! Define the phase space bounds (typical values for the Lorenz attractor)
  call geom%init(-30.0, 30.0, -30.0, 30.0, 0.0, 60.0)
  
  ! Print initial state
  print *, "Initial state:"
  call current_state%print()
  
  ! Check if initial state is within bounds
  if (geom%contains_point(current_state)) then
    print *, "Initial state is within phase space bounds."
  else
    print *, "Warning: Initial state is outside phase space bounds!"
  end if
  
  ! Integrate the system for some time steps
  num_steps = 1000
  do i = 1, num_steps
    call integrator%step(current_state)
    
    ! Print every 100 steps
    if (mod(i, 100) == 0) then
      print *, "Step ", i
      call current_state%print()
      
      ! Check if current state is within bounds
      if (.not. geom%contains_point(current_state)) then
        print *, "Warning: State has left the phase space bounds!"
      end if
    end if
  end do
  
  ! Compute distance from initial to final state using both methods
  ! Method 1: using the state's built-in distance method
  distance = initial_state%distance(current_state)
  print *, "Distance from initial to final state (state method): ", distance
  
  ! Method 2: using the module procedure for compatibility
  distance = compute_distance(initial_state, current_state)
  print *, "Distance from initial to final state (module procedure): ", distance
  
  ! Clean up resources
  call integrator%cleanup()
  
end program main 