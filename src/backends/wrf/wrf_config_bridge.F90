! Fortran bridge for WRF configuration system and domain management
! Provides C-callable wrappers for WRF initialization and domain lifecycle

module metada_wrf_config_bridge
  use iso_c_binding
  use module_configure, only: initial_config, model_to_grid_config_rec, &
                               grid_config_rec_type, model_config_rec, &
                               model_config_rec_type
  use module_domain, only: domain, head_grid, find_grid_by_id, &
                           alloc_and_configure_domain, dealloc_space_domain, &
                           init_module_domain
  implicit none

  ! Module-level storage for config_flags (not C-interoperable, so kept in Fortran)
  ! This is ALWAYS derived from WRF's authoritative model_config_rec via model_to_grid_config_rec
  ! NEVER manually update this variable - it must stay synchronized with WRF's configuration space
  ! TARGET attribute allows taking C pointers for passing to WRF observation operators
  type(grid_config_rec_type), save, target :: config_flags_
  
  ! Module-level flag to track WRF initialization state
  logical, save :: wrf_initialized_ = .false.

contains

  !============================================================================
  ! WRF Initialization Routines
  !============================================================================
  
  ! Initialize WRF modules (call once at startup)
  subroutine wrf_init_modules(phase) bind(C, name="wrf_init_modules_")
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
      ! Phase 1: Basic initialization
      call init_module_bc()
      call init_module_configure()
      call init_module_driver_constants()
      call init_module_model_constants()
      call init_module_domain()
      call init_module_machine()
      wrf_initialized_ = .true.
    else
      ! Phase 2: Advanced initialization (if needed)
      call init_module_wrf_error()
      call init_module_tiles()
    endif
    
  end subroutine wrf_init_modules
  
  ! Check if WRF has been initialized
  function wrf_is_initialized() bind(C, name="wrf_is_initialized_") result(initialized)
    implicit none
    logical(c_bool) :: initialized
    
    initialized = logical(wrf_initialized_, kind=c_bool)
  end function wrf_is_initialized

  ! C-callable wrapper for initial_config()
  ! Reads namelist.input and populates module-level model_config_rec
  subroutine wrf_initial_config() bind(C, name="wrf_initial_config_")
    implicit none
    
    ! Call WRF's initial_config to read namelist.input
    ! This populates the module-level model_config_rec structure
    call initial_config()
    
  end subroutine wrf_initial_config

  ! C-callable wrapper for model_to_grid_config_rec()
  ! Extracts domain-specific configuration from WRF's authoritative model_config_rec
  ! This is the ONLY way config_flags_ should be populated - ensuring we always
  ! use WRF's configuration space, not a parallel MetaDA configuration
  subroutine wrf_model_to_grid_config(domain_id) bind(C, name="wrf_model_to_grid_config_")
    implicit none
    integer(c_int), intent(in) :: domain_id
    
    ! Call WRF's model_to_grid_config_rec to extract domain-specific config
    ! This copies values from model_config_rec(domain_id) to module-level config_flags_
    ! All other code that needs to update config_flags_ must do so by:
    !   1. Updating model_config_rec first
    !   2. Then calling this subroutine to sync config_flags_
    call model_to_grid_config_rec(domain_id, model_config_rec, config_flags_)
    
  end subroutine wrf_model_to_grid_config

  !============================================================================
  ! WRF Domain Management Routines
  !============================================================================
  
  ! Allocate and initialize a WRF domain
  ! Returns 0 on success, non-zero on error
  function wrf_alloc_domain(domain_id) bind(C, name="wrf_alloc_domain_") result(ierr)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    integer(c_int) :: ierr
    type(domain), pointer :: parent
    logical :: active_this_task
    
    ierr = 0
    active_this_task = .true.
    nullify(parent)
    
    ! For domain 1 (root domain), parent is not used
    ! For nested domains, we would need to find the parent
    if (domain_id == 1) then
      ! Pass head_grid directly - alloc_and_configure_domain will allocate and set it
      call alloc_and_configure_domain(domain_id, active_this_task, head_grid, parent, -1)
    else
      ! For nested domains, find parent (parent_id = domain_id - 1 for simplicity)
      ! In a more complex setup, this would need proper parent tracking
      ierr = 1  ! Not yet implemented for nested domains
      return
    endif
    
    ! Verify that head_grid was successfully allocated
    if (.not. associated(head_grid)) then
      ierr = 1
    endif
    
  end function wrf_alloc_domain
  
  ! Deallocate a WRF domain
  subroutine wrf_dealloc_domain(domain_id) bind(C, name="wrf_dealloc_domain_")
    implicit none
    integer(c_int), intent(in), value :: domain_id
    
    call dealloc_space_domain(domain_id)
    
  end subroutine wrf_dealloc_domain
  
  ! Check if a domain exists
  function wrf_domain_exists(domain_id) bind(C, name="wrf_domain_exists_") result(exists)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    logical(c_bool) :: exists
    type(domain), pointer :: grid
    
    nullify(grid)
    call find_grid_by_id(domain_id, head_grid, grid)
    exists = logical(associated(grid), kind=c_bool)
    
  end function wrf_domain_exists
  
  !============================================================================
  ! WRF Grid Geometry Setters (write to WRF domain structures)
  !============================================================================
  
  ! Set grid geometry parameters in WRF structures
  ! Stores in WRF's model_config_rec (authoritative source) and grid domain structure
  ! Then syncs config_flags_ from model_config_rec to ensure consistency
  subroutine wrf_set_grid_geometry(domain_id, dx, dy, map_proj, &
                                    cen_lat, cen_lon, truelat1, truelat2, &
                                    stand_lon) bind(C, name="wrf_set_grid_geometry_")
    implicit none
    integer(c_int), intent(in), value :: domain_id
    real(c_float), intent(in), value :: dx, dy
    integer(c_int), intent(in), value :: map_proj
    real(c_float), intent(in), value :: cen_lat, cen_lon
    real(c_float), intent(in), value :: truelat1, truelat2, stand_lon
    
    type(domain), pointer :: grid
    
    ! 1. Store in WRF's model_config_rec (authoritative source for all WRF subsystems)
    !    This is indexed by domain_id and is the master configuration storage
    model_config_rec%dx(domain_id) = dx
    model_config_rec%dy(domain_id) = dy
    model_config_rec%map_proj(domain_id) = map_proj
    model_config_rec%cen_lat(domain_id) = cen_lat
    model_config_rec%cen_lon(domain_id) = cen_lon
    model_config_rec%truelat1(domain_id) = truelat1
    model_config_rec%truelat2(domain_id) = truelat2
    model_config_rec%stand_lon(domain_id) = stand_lon
    
    ! 2. Find and update the grid domain structure if it exists
    nullify(grid)
    call find_grid_by_id(domain_id, head_grid, grid)
    
    if (associated(grid)) then
      ! Set grid geometry in domain structure (for WRF routines that use grid%)
      grid%dx = dx
      grid%dy = dy
      grid%map_proj = map_proj
      grid%cen_lat = cen_lat
      grid%cen_lon = cen_lon
      grid%truelat1 = truelat1
      grid%truelat2 = truelat2
      grid%stand_lon = stand_lon
    endif
    
    ! 3. Sync config_flags_ from model_config_rec (NEVER manually update config_flags_)
    !    This ensures config_flags_ always reflects WRF's authoritative configuration
    call model_to_grid_config_rec(domain_id, model_config_rec, config_flags_)
    
  end subroutine wrf_set_grid_geometry
  
  !============================================================================
  ! WRF Grid Geometry Getters (read from WRF domain structures)
  !============================================================================
  
  ! Get grid geometry from WRF's model_config_rec (authoritative source)
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
  ! Grid Domain Pointer Access (for passing to WRF observation operators)
  !============================================================================
  
  ! Get C pointer to WRF grid/domain structure for a specific domain
  ! This allows METADA C++ code to pass grid to WRF Fortran routines
  function wrf_get_grid_ptr(domain_id) bind(C, name="wrf_get_grid_ptr_") result(ptr)
    implicit none
    integer(c_int), intent(in), value :: domain_id
    type(c_ptr) :: ptr
    type(domain), pointer :: grid
    
    nullify(grid)
    call find_grid_by_id(domain_id, head_grid, grid)
    
    if (associated(grid)) then
      ptr = c_loc(grid)
    else
      ptr = c_null_ptr
    endif
  end function wrf_get_grid_ptr
  
  !============================================================================
  ! Config Accessor Functions (for namelist values)
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

