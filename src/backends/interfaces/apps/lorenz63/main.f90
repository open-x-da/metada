!
! main.f90
! Main program to demonstrate the Lorenz 63 chaotic system
!

program main
  use state, only: state_type, state_init => init, state_print => print
  use geometry, only: geometry_type, geometry_init => init, compute_distance
  use model, only: model_type, model_init => init, step
  implicit none
  
  type(state_type) :: current_state, initial_state
  type(geometry_type) :: geom
  type(model_type) :: mdl
  
  integer :: i, num_steps
  real :: distance
  
  ! Initialize the system
  call state_init(current_state, 1.0, 1.0, 1.0)
  call state_init(initial_state, 1.0, 1.0, 1.0) ! Keep a copy of initial state
  
  ! Standard parameters for Lorenz 63 (sigma=10, rho=28, beta=8/3)
  call model_init(mdl, 10.0, 28.0, 8.0/3.0, 0.01)
  
  ! Define the phase space bounds (typical values for the Lorenz attractor)
  call geometry_init(geom, -30.0, 30.0, -30.0, 30.0, 0.0, 60.0)
  
  ! Print initial state
  print *, "Initial state:"
  call state_print(current_state)
  
  ! Integrate the system for some time steps
  num_steps = 1000
  do i = 1, num_steps
    call step(mdl, current_state)
    
    ! Print every 100 steps
    if (mod(i, 100) == 0) then
      print *, "Step ", i
      call state_print(current_state)
    end if
  end do
  
  ! Compute distance from initial to final state
  distance = compute_distance(initial_state, current_state)
  print *, "Distance from initial to final state: ", distance
  
end program main 