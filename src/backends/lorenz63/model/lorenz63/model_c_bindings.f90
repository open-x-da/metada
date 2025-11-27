!
! model_c_bindings.f90
! C bindings for Lorenz 63 model module
!

module model_c_bindings
  use, intrinsic :: iso_c_binding
  use model, only: lorenz63_model_type, rk4_integrator_type
  use state, only: state_type
  implicit none
  private
  
  ! Public interfaces for C bindings
  public :: c_lorenz63_model_create, c_lorenz63_model_destroy
  public :: c_lorenz63_model_get_params, c_lorenz63_model_set_params
  public :: c_rk4_integrator_create, c_rk4_integrator_destroy, c_rk4_integrator_step
  
contains

  !---------------------------------------------------------------------
  ! C-Fortran interoperability bindings
  !---------------------------------------------------------------------
  
  ! Create and return a C pointer to a new lorenz63 model
  function c_lorenz63_model_create(sigma, rho, beta, dt) bind(C, name="lorenz63_model_create") result(ptr_result)
    real(c_float), value, intent(in) :: sigma, rho, beta, dt
    type(c_ptr) :: ptr_result
    type(lorenz63_model_type), pointer :: f_ptr
    real :: sigma_f, rho_f, beta_f, dt_f
    
    allocate(f_ptr)
    ! Convert from C float to Fortran real (single precision)
    sigma_f = real(sigma)
    rho_f = real(rho)
    beta_f = real(beta)
    dt_f = real(dt)
    call f_ptr%init(sigma_f, rho_f, beta_f, dt_f)
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
    real :: sigma_f, rho_f, beta_f, dt_f
    
    call c_f_pointer(model_ptr, f_ptr)
    sigma_f = f_ptr%get_sigma()
    rho_f = f_ptr%get_rho()
    beta_f = f_ptr%get_beta()
    dt_f = f_ptr%get_time_step()
    ! Convert from Fortran real to C float (both single precision)
    sigma = real(sigma_f, c_float)
    rho = real(rho_f, c_float)
    beta = real(beta_f, c_float)
    dt = real(dt_f, c_float)
  end subroutine c_lorenz63_model_get_params
  
  ! Set the parameters of a lorenz63 model
  subroutine c_lorenz63_model_set_params(model_ptr, sigma, rho, beta, dt) bind(C, name="lorenz63_model_set_params")
    type(c_ptr), value, intent(in) :: model_ptr
    real(c_float), value, intent(in) :: sigma, rho, beta, dt
    type(lorenz63_model_type), pointer :: f_ptr
    real :: sigma_f, rho_f, beta_f, dt_f
    
    call c_f_pointer(model_ptr, f_ptr)
    ! Convert from C float to Fortran real (single precision)
    sigma_f = real(sigma)
    rho_f = real(rho)
    beta_f = real(beta)
    dt_f = real(dt)
    call f_ptr%init(sigma_f, rho_f, beta_f, dt_f)
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
  
end module model_c_bindings 