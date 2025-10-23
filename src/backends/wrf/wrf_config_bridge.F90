! Fortran bridge for WRFDA configuration system and domain management
! Provides C-callable wrappers for WRFDA initialization following the standard workflow
! from da_wrfvar_init1.inc and da_wrfvar_init2.inc

module metada_wrf_config_bridge
  use iso_c_binding
  use module_configure, only: initial_config, model_to_grid_config_rec, &
                               grid_config_rec_type, model_config_rec, &
                               model_config_rec_type
  use module_domain, only: domain, head_grid, alloc_and_configure_domain
  use da_wrf_interfaces, only: set_scalar_indices_from_config, setup_timekeeping, &
                                init_wrfio
  implicit none

  ! Module-level storage for config_flags (not C-interoperable, so kept in Fortran)
  ! This is ALWAYS derived from WRFDA's authoritative model_config_rec via model_to_grid_config_rec
  ! NEVER manually update this variable - it must stay synchronized with WRFDA's configuration space
  ! TARGET attribute allows taking C pointers for passing to WRFDA observation operators
  type(grid_config_rec_type), save, target :: config_flags_
  
  ! Module-level pointer to null domain for grid allocation
  type(domain), pointer, save :: null_domain => null()
  
  ! Module-level flag to track WRFDA initialization state
  logical, save :: wrfda_initialized_ = .false.

contains

  !============================================================================
  ! WRFDA Initialization Routines
  !============================================================================
  
  ! Initialize WRFU (WRF ESMF time utilities)
  subroutine wrfda_wrfu_initialize() bind(C, name="wrfda_wrfu_initialize_")
    use module_utility, only: WRFU_Initialize, WRFU_CAL_GREGORIAN
    implicit none
    
    ! Initialize WRFU time manager with Gregorian calendar
    ! Note: Could add support for NO_LEAP_CALENDAR via a parameter if needed
    call WRFU_Initialize(defaultCalKind=WRFU_CAL_GREGORIAN)
    
  end subroutine wrfda_wrfu_initialize
  
  ! Initialize WRFDA modules (call once at startup)
  ! Note: Domain management is handled separately by WRFDA itself
  subroutine wrfda_init_modules(phase) bind(C, name="wrfda_init_modules_")
    use module_bc, only: init_module_bc
    use module_configure, only: init_module_configure
    use module_driver_constants, only: init_module_driver_constants
    use module_model_constants, only: init_module_model_constants
    use module_machine, only: init_module_machine
    use module_tiles, only: init_module_tiles
    use module_wrf_error, only: init_module_wrf_error
    implicit none
    integer(c_int), intent(in), value :: phase
    
    if (phase == 1) then
      ! Phase 1: Basic initialization (configuration and constants only)
      call init_module_bc()
      call init_module_configure()
      call init_module_driver_constants()
      call init_module_model_constants()
      call init_module_machine()
      wrfda_initialized_ = .true.
    else
      ! Phase 2: Advanced initialization (if needed)
      call init_module_wrf_error()
      call init_module_tiles()
    endif
    
  end subroutine wrfda_init_modules
  
  ! Check if WRFDA has been initialized
  function wrfda_is_initialized() bind(C, name="wrfda_is_initialized_") result(initialized)
    implicit none
    logical(c_bool) :: initialized
    
    initialized = logical(wrfda_initialized_, kind=c_bool)
  end function wrfda_is_initialized

  ! C-callable wrapper for initial_config()
  ! Reads namelist.input and populates module-level model_config_rec
  subroutine wrf_initial_config() bind(C, name="wrf_initial_config_")
    implicit none
    
    ! Call WRF's initial_config to read namelist.input
    ! This populates the module-level model_config_rec structure
    call initial_config()
    
  end subroutine wrf_initial_config

  ! C-callable wrapper for model_to_grid_config_rec()
  ! Extracts domain-specific configuration from WRFDA's authoritative model_config_rec
  ! This is the ONLY way config_flags_ should be populated - ensuring we always
  ! use WRFDA's configuration space, not a parallel MetaDA configuration
  subroutine wrf_model_to_grid_config(domain_id) bind(C, name="wrf_model_to_grid_config_")
    implicit none
    integer(c_int), intent(in) :: domain_id
    
    ! Call WRFDA's model_to_grid_config_rec to extract domain-specific config
    ! This copies values from model_config_rec(domain_id) to module-level config_flags_
    ! All other code that needs to update config_flags_ must do so by:
    !   1. Updating model_config_rec first
    !   2. Then calling this subroutine to sync config_flags_
    call model_to_grid_config_rec(domain_id, model_config_rec, config_flags_)
    
  end subroutine wrf_model_to_grid_config

  !============================================================================
  ! WRFDA Domain Allocation (following da_wrfvar_init2.inc workflow)
  !============================================================================
  
  !> @brief Allocate and configure WRFDA domain following standard workflow
  !> @details Follows the exact sequence from da_wrfvar_init2.inc:
  !>   1. alloc_and_configure_domain (allocates head_grid)
  !>   2. model_to_grid_config_rec (extracts domain config)
  !>   3. set_scalar_indices_from_config (sets tracer indices)
  !>   4. init_wrfio (initializes WRF I/O system)
  !>   5. setup_timekeeping (sets up domain clock)
  !> @param[in] domain_id Domain ID (typically 1 for single-domain)
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_alloc_and_init_domain(domain_id) &
      bind(C, name="wrfda_alloc_and_init_domain_") result(ierr)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    integer :: idum1, idum2
    
    ierr = 0
    
    ! Ensure null_domain is nullified
    if (.not. associated(null_domain)) then
      nullify(null_domain)
    end if
    
    ! STEP 1: Allocate and configure domain (WRFDA standard workflow)
    call alloc_and_configure_domain(domain_id=domain_id, &
                                     grid=head_grid, &
                                     parent=null_domain, &
                                     kid=-1)
    
    ! Verify allocation succeeded
    if (.not. associated(head_grid)) then
      ierr = 1
      return
    end if
    
    ! STEP 2: Extract domain-specific config
    call model_to_grid_config_rec(head_grid%id, model_config_rec, config_flags_)
    
    ! STEP 3: Set scalar indices from config
    call set_scalar_indices_from_config(head_grid%id, idum1, idum2)
    
    ! STEP 4: Initialize WRF I/O system
    call init_wrfio()
    
    ! STEP 5: Setup timekeeping for the domain
    call setup_timekeeping(head_grid)
    
  end function wrfda_alloc_and_init_domain
  
  !> @brief Get pointer to WRFDA's head_grid
  !> @return C pointer to head_grid, or c_null_ptr if not allocated
  type(c_ptr) function wrfda_get_head_grid_ptr() bind(C, name="wrfda_get_head_grid_ptr_")
    implicit none
    
    if (associated(head_grid)) then
      wrfda_get_head_grid_ptr = c_loc(head_grid)
    else
      wrfda_get_head_grid_ptr = c_null_ptr
    end if
    
  end function wrfda_get_head_grid_ptr
  
  !> @brief Check if head_grid has been allocated
  !> @return True if head_grid is allocated
  logical(c_bool) function wrfda_head_grid_allocated() bind(C, name="wrfda_head_grid_allocated_")
    implicit none
    
    wrfda_head_grid_allocated = logical(associated(head_grid), kind=c_bool)
    
  end function wrfda_head_grid_allocated

  !============================================================================
  ! WRFDA Grid Geometry Getters (read from model_config_rec)
  !============================================================================
  
  ! Get grid geometry from WRFDA's model_config_rec (authoritative source)
  ! These getters read from model_config_rec which is the master configuration storage
  
  function wrf_get_grid_dx(domain_id) bind(C, name="wrf_get_grid_dx_") result(dx)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    real(c_float) :: dx
    
    dx = real(model_config_rec%dx(domain_id), kind=c_float)
  end function wrf_get_grid_dx
  
  function wrf_get_grid_dy(domain_id) bind(C, name="wrf_get_grid_dy_") result(dy)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    real(c_float) :: dy
    
    dy = real(model_config_rec%dy(domain_id), kind=c_float)
  end function wrf_get_grid_dy
  
  function wrf_get_grid_map_proj(domain_id) bind(C, name="wrf_get_grid_map_proj_") result(map_proj)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    integer(c_int) :: map_proj
    
    map_proj = int(model_config_rec%map_proj(domain_id), kind=c_int)
  end function wrf_get_grid_map_proj
  
  function wrf_get_grid_cen_lat(domain_id) bind(C, name="wrf_get_grid_cen_lat_") result(cen_lat)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    real(c_float) :: cen_lat
    
    cen_lat = real(model_config_rec%cen_lat(domain_id), kind=c_float)
  end function wrf_get_grid_cen_lat
  
  function wrf_get_grid_cen_lon(domain_id) bind(C, name="wrf_get_grid_cen_lon_") result(cen_lon)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    real(c_float) :: cen_lon
    
    cen_lon = real(model_config_rec%cen_lon(domain_id), kind=c_float)
  end function wrf_get_grid_cen_lon
  
  function wrf_get_grid_truelat1(domain_id) bind(C, name="wrf_get_grid_truelat1_") result(truelat1)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    real(c_float) :: truelat1
    
    truelat1 = real(model_config_rec%truelat1(domain_id), kind=c_float)
  end function wrf_get_grid_truelat1
  
  function wrf_get_grid_truelat2(domain_id) bind(C, name="wrf_get_grid_truelat2_") result(truelat2)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    real(c_float) :: truelat2
    
    truelat2 = real(model_config_rec%truelat2(domain_id), kind=c_float)
  end function wrf_get_grid_truelat2
  
  function wrf_get_grid_stand_lon(domain_id) bind(C, name="wrf_get_grid_stand_lon_") result(stand_lon)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    real(c_float) :: stand_lon
    
    stand_lon = real(model_config_rec%stand_lon(domain_id), kind=c_float)
  end function wrf_get_grid_stand_lon
  
  !============================================================================
  ! Config Flags Pointer Access (for passing to WRF observation operators)
  !============================================================================
  
  ! Get C pointer to config_flags_ for passing to WRF observation operators
  ! This allows METADA C++ code to pass config_flags to WRF Fortran routines
  function wrf_get_config_flags_ptr() bind(C, name="wrf_get_config_flags_ptr_") result(ptr)
    implicit none
    type(c_ptr) :: ptr
    
    ! Return C pointer to module-level config_flags_
    ptr = c_loc(config_flags_)
  end function wrf_get_config_flags_ptr
  
  ! Get size of config_flags structure (for validation)
  function wrf_get_config_flags_size() bind(C, name="wrf_get_config_flags_size_") result(size_bytes)
    implicit none
    integer(c_size_t) :: size_bytes
    
    ! Calculate size of config_flags_ structure
    size_bytes = storage_size(config_flags_) / 8  ! Convert bits to bytes
  end function wrf_get_config_flags_size
  
  !============================================================================
  ! Config Accessor Functions (for namelist values from config_flags_)
  !============================================================================
  ! Accessor functions for config_flags_ members
  ! These access the module-level config_flags_ variable
  
  function wrf_get_dx() bind(C, name="wrf_get_dx_") result(dx)
    implicit none
    real(c_double) :: dx
    
    dx = real(config_flags_%dx, kind=c_double)
  end function wrf_get_dx

  function wrf_get_dy() bind(C, name="wrf_get_dy_") result(dy)
    implicit none
    real(c_double) :: dy
    
    dy = real(config_flags_%dy, kind=c_double)
  end function wrf_get_dy

  function wrf_get_map_proj() bind(C, name="wrf_get_map_proj_") result(map_proj)
    implicit none
    integer(c_int) :: map_proj
    
    map_proj = int(config_flags_%map_proj, kind=c_int)
  end function wrf_get_map_proj

  function wrf_get_truelat1() bind(C, name="wrf_get_truelat1_") result(truelat1)
    implicit none
    real(c_double) :: truelat1
    
    truelat1 = real(config_flags_%truelat1, kind=c_double)
  end function wrf_get_truelat1

  function wrf_get_truelat2() bind(C, name="wrf_get_truelat2_") result(truelat2)
    implicit none
    real(c_double) :: truelat2
    
    truelat2 = real(config_flags_%truelat2, kind=c_double)
  end function wrf_get_truelat2

  function wrf_get_stand_lon() bind(C, name="wrf_get_stand_lon_") result(stand_lon)
    implicit none
    real(c_double) :: stand_lon
    
    stand_lon = real(config_flags_%stand_lon, kind=c_double)
  end function wrf_get_stand_lon

  function wrf_get_cen_lat() bind(C, name="wrf_get_cen_lat_") result(cen_lat)
    implicit none
    real(c_double) :: cen_lat
    
    cen_lat = real(config_flags_%cen_lat, kind=c_double)
  end function wrf_get_cen_lat

  function wrf_get_cen_lon() bind(C, name="wrf_get_cen_lon_") result(cen_lon)
    implicit none
    real(c_double) :: cen_lon
    
    cen_lon = real(config_flags_%cen_lon, kind=c_double)
  end function wrf_get_cen_lon

end module metada_wrf_config_bridge

