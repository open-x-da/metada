! Dispatch C-bindings that translate flat arrays into WRFDA types
! and call the included transform routines.
!
! OUT-OF-DOMAIN OBSERVATION HANDLING:
! This module implements WRFDA-compliant handling of out-of-domain observations:
! 1. find_fractional_ij returns -1 for i,j and -1.0 for xfloat,yfloat when observations are out of domain
! 2. Callers check for -1 values and handle them gracefully by:
!    - Setting invalid indices (-1) in the iv structure
!    - Skipping further processing for that observation
!    - Continuing with the next observation
! 3. This follows WRFDA's standard pattern of setting outside=.true. and returning immediately

module metada_wrfda_dispatch
  use iso_c_binding
  use module_domain,        only: domain, x_type
  use da_define_structures, only: iv_type, y_type, xbx_type, da_allocate_y, da_allocate_obs_info, da_zero_y
  use module_symbols_util,  only: wrfu_initialize, wrfu_finalize, wrfu_cal_gregorian
  use module_symbols_util, only: WRFU_ClockCreate, WRFU_TimeIntervalSet, WRFU_SUCCESS, WRFU_Time, WRFU_TimeInterval, WRFU_INITIALIZE, WRFU_CAL_GREGORIAN
  use da_control, only: sound, synop, pilot, satem, geoamv, polaramv, airep, gpspw, gpsref, &
                        metar, ships, ssmi_rv, ssmi_tb, ssmt1, ssmt2, qscat, profiler, buoy, &
                        bogus, pseudo, radar, radiance, airsr, sonde_sfc, mtgirs, tamdar, &
                        tamdar_sfc, rain, gpseph, lightning, &
                        sfc_assi_options, sfc_assi_options_1, trace_use_dull, &
                        var4d_run, num_fgat_time, missing_r, missing_data, num_ob_indexes, &
                        kts, kte, its, ite, jts, jte, Max_StHeight_Diff, kms, kme, &
                        ids, ide, jds, jde, ims, ime, jms, jme, ips, ipe, jps, jpe, kds, kde, kps, kpe, &
                        myproc, num_procs, rootproc, comm, num_qcstat_conv
  use da_tools, only: proj_info, da_map_set, da_llxy_wrf, da_togrid
  use da_metar,  only: da_transform_xtoy_metar,  da_transform_xtoy_metar_adj
  use da_synop,  only: da_transform_xtoy_synop,  da_transform_xtoy_synop_adj
  use da_buoy,   only: da_transform_xtoy_buoy,   da_transform_xtoy_buoy_adj
  use da_ships,  only: da_transform_xtoy_ships,  da_transform_xtoy_ships_adj
  use da_airep,  only: da_transform_xtoy_airep,  da_transform_xtoy_airep_adj
  use da_pilot,  only: da_transform_xtoy_pilot,  da_transform_xtoy_pilot_adj
  use da_sound,  only: da_transform_xtoy_sound,  da_transform_xtoy_sound_adj, &
                       da_transform_xtoy_sonde_sfc, da_transform_xtoy_sonde_sfc_adj
  use da_par_util, only: da_copy_dims, da_copy_tile_dims
  use da_tools, only: da_togrid, da_llxy_wrf, da_map_set, da_map_init
  use module_configure, only: grid_config_rec_type
  use da_minimisation, only: da_get_innov_vector

  implicit none
  
  ! Grid pointer: points to grid allocated via wrfda_alloc_and_init_domain
  ! This allows observation operators to access the grid
  type(domain), pointer, save :: persistent_grid => null()
  logical, save :: grid_initialized = .false.
  
  ! Module-level observation structures (allocated by wrfda_read_and_allocate_observations)
  ! Note: With intent(inout) in da_setup_obs_structures, module-level variables work correctly
  type(iv_type), save, target :: wrfda_iv
  type(y_type), save, target :: wrfda_ob
  logical, save :: wrfda_obs_allocated = .false.
  
  ! Temporary y_type for tangent linear output (to avoid corrupting wrfda_ob)
  ! Allocated on first TL call, deallocated on cleanup
  type(y_type), pointer, save :: wrfda_y_tl => null()
  logical, save :: wrfda_y_tl_allocated = .false.
  
  ! Map projection information for WRFDA coordinate conversion
  ! Define projection constants
  integer, parameter :: PROJ_LATLON = 0
  integer, parameter :: PROJ_MERC = 1
  integer, parameter :: PROJ_PS = 2
  integer, parameter :: PROJ_LC = 3
  
  ! Define proj_info type (simplified version of WRFDA's proj_info)
  ! Use WRFDA's proj_info type from da_tools module
  
  ! Module-level map_info structure initialized by da_setup_firstguess_wrf → da_map_set
  type(proj_info), save :: map_info
  
  ! Persistent y_type for weighted residual (jo_grad_y = -R^{-1} · residual)
  ! Computed by wrfda_compute_weighted_residual, used by wrfda_xtoy_adjoint_grid
  type(y_type), save :: persistent_jo_grad_y
  
contains


  !> @brief Tangent linear operator: H'(xb)·δx using WRFDA's proven da_transform_xtoy
  !> @details Uses WRFDA's top-level da_transform_xtoy which automatically dispatches
  !>          to all observation types (sound, synop, metar, ships, radar, satellite, etc.)
  !>          based on iv%info(*)%nlocal counts.
  !>
  !>          In WRFDA's incremental formulation:
  !>          - xb (grid%xb): background state, fixed during inner loop
  !>          - δx (grid%xa): analysis increment, updated by minimizer
  !>          - This computes H'(xb)·δx → y using WRFDA's proven dispatch logic
  !>
  !> @note Output extraction is delegated to C++ code which calls
  !>       wrfda_extract_observations() for consistency with nonlinear operator.
  integer(c_int) function wrfda_xtoy_apply_grid() bind(C, name="wrfda_xtoy_apply_grid")
    use module_domain, only: head_grid
    use da_obs, only: da_transform_xtoy
    use da_define_structures, only: da_allocate_y
    implicit none
    
    real :: dummy_cv(1)  ! Dummy control variable (not used in TL operator)
    
    ! Validate module-level structures
    if (.not. associated(head_grid)) then
      wrfda_xtoy_apply_grid = 1_c_int
      return
    end if
    
    if (.not. wrfda_obs_allocated) then
      wrfda_xtoy_apply_grid = 1_c_int
      return
    end if
    
    ! Allocate temporary y_type for TL output (preserves wrfda_ob)
    ! Only allocate once, reuse for subsequent TL calls
    if (.not. wrfda_y_tl_allocated) then
      allocate(wrfda_y_tl)
      call da_allocate_y(wrfda_iv, wrfda_y_tl)
      wrfda_y_tl_allocated = .true.
    end if
    
    ! Use WRFDA's proven top-level function which automatically dispatches
    ! to ALL observation types (sound, synop, metar, ships, buoy, pilot, airep,
    ! gpsref, radar, satellite, etc.) based on iv%info(*)%nlocal
    ! CRITICAL: Use wrfda_y_tl instead of wrfda_ob to preserve observation values
    ! Note: cv_size and cv are not needed for TL operator (only for full transform)
    call da_transform_xtoy(0, dummy_cv, head_grid, wrfda_iv, wrfda_y_tl)
    
    wrfda_xtoy_apply_grid = 0_c_int
  end function wrfda_xtoy_apply_grid

  !> @brief Compute weighted residual using WRFDA's proven workflow
  !> @details Replaces MetaDA's manual residual and R^{-1} computation with WRFDA's
  !>          proven functions:
  !>          1. da_calculate_residual(iv, y, re): re = (O-B) - H(xa)
  !>          2. da_calculate_grady(iv, re, jo_grad_y): jo_grad_y = -R^{-1} · re
  !>
  !>          The result is stored in module-level persistent_delta_y for use by
  !>          the adjoint operator.
  !>
  !> @param[in] iv_ptr Innovation vector (O-B with error statistics)
  !> @param[in] y_ptr Simulated observations H(xa) or H'(δx)
  !> @return 0 on success, non-zero on error
  integer(c_int) function wrfda_compute_weighted_residual(iv_ptr, y_ptr) &
      bind(C, name="wrfda_compute_weighted_residual")
    use da_minimisation, only: da_calculate_residual, da_calculate_grady
    use da_define_structures, only: da_allocate_y, da_deallocate_y
    implicit none
    type(c_ptr), value :: iv_ptr, y_ptr
    
    type(iv_type), pointer :: iv
    type(y_type), pointer :: y  
    type(y_type), target :: re, jo_grad_y
    
    call c_f_pointer(iv_ptr, iv)
    call c_f_pointer(y_ptr, y)
    
    ! Allocate residual and gradient structures
    call da_allocate_y(iv, re)
    call da_allocate_y(iv, jo_grad_y)
    
    ! Step 1: WRFDA's proven residual calculation
    ! re = (O-B) - H(xa) for ALL observation types
    call da_calculate_residual(iv, y, re)
    
    ! Step 2: WRFDA's proven gradient calculation
    ! jo_grad_y = -R^{-1} · re for ALL observation types
    call da_calculate_grady(iv, re, jo_grad_y)
    
    ! Step 3: Store jo_grad_y in persistent module variable
    ! Deallocate previous if exists
    if (associated(persistent_jo_grad_y%synop)) call da_deallocate_y(persistent_jo_grad_y)
    if (associated(persistent_jo_grad_y%sound)) call da_deallocate_y(persistent_jo_grad_y)
    
    ! Transfer ownership to persistent variable (no copy needed)
    persistent_jo_grad_y = jo_grad_y
    
    ! Cleanup residual only (jo_grad_y is now in persistent storage)
    call da_deallocate_y(re)
    
    wrfda_compute_weighted_residual = 0_c_int
  end function wrfda_compute_weighted_residual

  !> @brief Adjoint operator: H^T·jo_grad_y using WRFDA's proven da_transform_xtoy_adj
  !> @details Uses WRFDA's top-level da_transform_xtoy_adj which automatically dispatches
  !>          to all observation types based on iv%info(*)%nlocal counts.
  !>
  !>          Computes H^T·jo_grad_y → grid%xa where:
  !>          - jo_grad_y: weighted residual from wrfda_compute_weighted_residual
  !>          - grid%xa: state space gradient (output)
  !>
  !> @note Uses persistent_jo_grad_y computed by wrfda_compute_weighted_residual.
  !>       Output is written directly to grid%xa by WRFDA.
  integer(c_int) function wrfda_xtoy_adjoint_grid(grid_ptr, iv_ptr) bind(C, name="wrfda_xtoy_adjoint_grid")
    use da_obs, only: da_transform_xtoy_adj
    use da_define_structures, only: da_zero_x
    implicit none
    type(c_ptr), value :: grid_ptr, iv_ptr

    type(domain), pointer :: grid
    type(iv_type), pointer :: iv
    real :: dummy_cv(1)  ! Dummy control variable (not used in adjoint operator)
    
    ! Convert C pointers to Fortran pointers
    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(iv_ptr, iv)
    
    if (.not. associated(iv)) then
      wrfda_xtoy_adjoint_grid = 1_c_int
      return
    end if
    
    ! Zero ALL analysis increment fields using WRFDA's proven function
    call da_zero_x(grid%xa)
    
    ! Use WRFDA's proven top-level adjoint function with persistent jo_grad_y
    ! persistent_jo_grad_y was computed by wrfda_compute_weighted_residual using:
    !   1. da_calculate_residual: re = (O-B) - H(xa)
    !   2. da_calculate_grady: jo_grad_y = -R^{-1} · re
    ! This automatically handles ALL observation types (sound, synop, metar, etc.)
    call da_transform_xtoy_adj(0, dummy_cv, grid, iv, persistent_jo_grad_y, grid%xa)
    
    ! grid%xa now contains H^T · jo_grad_y for ALL observation types!
    
    wrfda_xtoy_adjoint_grid = 0_c_int
  end function wrfda_xtoy_adjoint_grid

  ! Copy gradient from jo_grad_x to state arrays (gradient accumulation)
  ! This function accumulates the adjoint gradients into the state space arrays
  ! For incremental 3D-Var, this represents the gradient of the cost function with respect to state
  subroutine copy_x_to_state(jo_grad_x, u, v, t, q, psfc, nx, ny, nz)
    type(x_type), intent(in) :: jo_grad_x
    real(c_double), intent(inout) :: u(*), v(*), t(*), q(*), psfc(*)
    integer, intent(in) :: nx, ny, nz
    integer :: i,j,k, nz1
    nz1 = max(1, nz)
    
    ! Accumulate gradients into state arrays (gradient accumulation)
    do k=1,nz1; do j=1,ny; do i=1,nx
      u(i + (j-1)*nx + (k-1)*nx*ny) = u(i + (j-1)*nx + (k-1)*nx*ny) + real(jo_grad_x%u(i,j,k), kind=c_double)
      v(i + (j-1)*nx + (k-1)*nx*ny) = v(i + (j-1)*nx + (k-1)*nx*ny) + real(jo_grad_x%v(i,j,k), kind=c_double)
      t(i + (j-1)*nx + (k-1)*nx*ny) = t(i + (j-1)*nx + (k-1)*nx*ny) + real(jo_grad_x%t(i,j,k), kind=c_double)
      q(i + (j-1)*nx + (k-1)*nx*ny) = q(i + (j-1)*nx + (k-1)*nx*ny) + real(jo_grad_x%q(i,j,k), kind=c_double)
    end do; end do; end do
    do j=1,ny; do i=1,nx
      psfc(i + (j-1)*nx) = psfc(i + (j-1)*nx) + real(jo_grad_x%psfc(i,j), kind=c_double)
    end do; end do
  end subroutine copy_x_to_state

  subroutine zero_x_like(x, nx, ny, nz)
    type(x_type), intent(inout) :: x
    integer, intent(in) :: nx, ny, nz
    integer :: nz1
    nz1 = max(1, nz)
    allocate(x%u(nx,ny,nz1)); x%u = 0.0
    allocate(x%v(nx,ny,nz1)); x%v = 0.0
    allocate(x%t(nx,ny,nz1)); x%t = 0.0
    allocate(x%q(nx,ny,nz1)); x%q = 0.0
    allocate(x%psfc(nx,ny));  x%psfc = 0.0
  end subroutine zero_x_like

  !============================================================================
  ! WRFDA Native Increment Operations (using proven WRFDA functions)
  !============================================================================
  
  ! Copy analysis increments from one grid to another using WRFDA's da_copy_xa
  ! This copies ALL fields in x_type (35+ fields), not just control variables
  subroutine wrfda_copy_xa(grid_dst_ptr, grid_src_ptr) bind(C, name="wrfda_copy_xa")
    use da_vtox_transforms, only : da_copy_xa
    implicit none
    type(c_ptr), value :: grid_dst_ptr, grid_src_ptr
    type(domain), pointer :: grid_dst, grid_src
    
    call c_f_pointer(grid_dst_ptr, grid_dst)
    call c_f_pointer(grid_src_ptr, grid_src)
    
    ! Use WRFDA's proven function to copy ALL increment fields
    call da_copy_xa(grid_dst%xa, grid_src%xa)
  end subroutine wrfda_copy_xa
  
  ! Randomize ALL analysis increment fields (following da_copy_xa structure)
  ! This ensures all 35+ fields are randomized, not just control variables
  subroutine wrfda_randomize_xa(grid_ptr) bind(C, name="wrfda_randomize_xa")
    use da_define_structures, only : da_zero_x
    implicit none
    type(c_ptr), value :: grid_ptr
    type(domain), pointer :: grid
    
    call c_f_pointer(grid_ptr, grid)
    
    ! First zero all fields using WRFDA's proven function
    call da_zero_x(grid%xa)
    
    ! Randomize ALL fields using array intrinsics (following da_copy_xa structure)
    ! This uses the actual allocated bounds of the arrays
    
    ! 3D fields - Primary variables
    call random_number(grid%xa%u); grid%xa%u = grid%xa%u - 0.5
    call random_number(grid%xa%v); grid%xa%v = grid%xa%v - 0.5
    call random_number(grid%xa%t); grid%xa%t = grid%xa%t - 0.5
    call random_number(grid%xa%q); grid%xa%q = grid%xa%q - 0.5
    
    ! 3D fields - Hydrometeors
    call random_number(grid%xa%qcw); grid%xa%qcw = grid%xa%qcw - 0.5
    call random_number(grid%xa%qrn); grid%xa%qrn = grid%xa%qrn - 0.5
    call random_number(grid%xa%qci); grid%xa%qci = grid%xa%qci - 0.5
    call random_number(grid%xa%qsn); grid%xa%qsn = grid%xa%qsn - 0.5
    call random_number(grid%xa%qgr); grid%xa%qgr = grid%xa%qgr - 0.5
    
    ! 3D fields - Diagnostic variables
    call random_number(grid%xa%w); grid%xa%w = grid%xa%w - 0.5
    call random_number(grid%xa%p); grid%xa%p = grid%xa%p - 0.5
    call random_number(grid%xa%geoh); grid%xa%geoh = grid%xa%geoh - 0.5
    call random_number(grid%xa%rh); grid%xa%rh = grid%xa%rh - 0.5
    call random_number(grid%xa%wh); grid%xa%wh = grid%xa%wh - 0.5
    call random_number(grid%xa%rho); grid%xa%rho = grid%xa%rho - 0.5
    call random_number(grid%xa%ref); grid%xa%ref = grid%xa%ref - 0.5
    call random_number(grid%xa%qt); grid%xa%qt = grid%xa%qt - 0.5
    
    ! 2D fields - Surface variables
    call random_number(grid%xa%psfc); grid%xa%psfc = grid%xa%psfc - 0.5
    call random_number(grid%xa%mu); grid%xa%mu = grid%xa%mu - 0.5
    call random_number(grid%xa%tgrn); grid%xa%tgrn = grid%xa%tgrn - 0.5
    call random_number(grid%xa%u10); grid%xa%u10 = grid%xa%u10 - 0.5
    call random_number(grid%xa%v10); grid%xa%v10 = grid%xa%v10 - 0.5
    call random_number(grid%xa%t2); grid%xa%t2 = grid%xa%t2 - 0.5
    call random_number(grid%xa%q2); grid%xa%q2 = grid%xa%q2 - 0.5
    
    ! 2D fields - Derived variables
    call random_number(grid%xa%ztd); grid%xa%ztd = grid%xa%ztd - 0.5
    call random_number(grid%xa%tpw); grid%xa%tpw = grid%xa%tpw - 0.5
    call random_number(grid%xa%speed); grid%xa%speed = grid%xa%speed - 0.5
    
    ! 2D fields - Brightness temperatures
    call random_number(grid%xa%tb19v); grid%xa%tb19v = grid%xa%tb19v - 0.5
    call random_number(grid%xa%tb19h); grid%xa%tb19h = grid%xa%tb19h - 0.5
    call random_number(grid%xa%tb22v); grid%xa%tb22v = grid%xa%tb22v - 0.5
    call random_number(grid%xa%tb37v); grid%xa%tb37v = grid%xa%tb37v - 0.5
    call random_number(grid%xa%tb37h); grid%xa%tb37h = grid%xa%tb37h - 0.5
    call random_number(grid%xa%tb85v); grid%xa%tb85v = grid%xa%tb85v - 0.5
    call random_number(grid%xa%tb85h); grid%xa%tb85h = grid%xa%tb85h - 0.5
  end subroutine wrfda_randomize_xa

  !============================================================================
  ! 5-variable interface (control variables only)
  !============================================================================
  
  ! Update analysis increments (xa) using WRFDA's proven da_zero_x function
  ! This ensures ALL increment fields are properly initialized, not just the 5 control variables
  ! NOTE: For complete x_type alignment, use wrfda_copy_xa instead
  !> @brief Update all analysis increment fields in grid%xa from C++ arrays
  !> @details Updates all 18 identified fields in grid%xa:
  !>          - 3D fields: u, v, w, t, p, q, qt, rh, rho, geoh, wh, qcw, qrn, qci, qsn, qgr
  !>          - 2D fields: psfc, mu
  !> @param[in] u_inc, v_inc, w_inc, t_inc, p_inc, q_inc, qt_inc, rh_inc, rho_inc, geoh_inc, wh_inc
  !>            qcw_inc, qrn_inc, qci_inc, qsn_inc, qgr_inc 3D arrays (size: nx*ny*nz)
  !> @param[in] psfc_inc, mu_inc 2D arrays (size: nx*ny)
  !> @param[in] grid_ptr Pointer to WRFDA domain structure
  subroutine wrfda_update_analysis_increments(u_inc, v_inc, w_inc, t_inc, p_inc, q_inc, &
                                              qt_inc, rh_inc, rho_inc, geoh_inc, wh_inc, &
                                              qcw_inc, qrn_inc, qci_inc, qsn_inc, qgr_inc, &
                                              psfc_inc, mu_inc, grid_ptr) &
      bind(C, name="wrfda_update_analysis_increments")
    use da_define_structures, only : da_zero_x
    implicit none
    ! 3D field arrays
    real(c_double), intent(in) :: u_inc(*), v_inc(*), w_inc(*), t_inc(*), p_inc(*), q_inc(*)
    real(c_double), intent(in) :: qt_inc(*), rh_inc(*), rho_inc(*), geoh_inc(*), wh_inc(*)
    real(c_double), intent(in) :: qcw_inc(*), qrn_inc(*), qci_inc(*), qsn_inc(*), qgr_inc(*)
    ! 2D field arrays
    real(c_double), intent(in) :: psfc_inc(*), mu_inc(*)
    type(c_ptr), value :: grid_ptr
    
    type(domain), pointer :: grid
    integer :: i, j, k, nx, ny, nz, idx
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(grid_ptr, grid)
    
    ! Use WRFDA's proven function to zero ALL increment fields
    ! This properly initializes all 50+ fields in x_type
    call da_zero_x(grid%xa)
    
    ! Get grid dimensions
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    
    ! Populate 3D fields from C++ arrays (row-major indexing)
    ! Note: Some fields may be dummy arrays (1x1x1) if not used (e.g., cloud_cv_options <= 1)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          idx = i + (j-1)*nx + (k-1)*nx*ny
          ! Primary state variables (always allocated with full dimensions)
          if (associated(grid%xa%u)) grid%xa%u(i,j,k) = real(u_inc(idx), kind=4)
          if (associated(grid%xa%v)) grid%xa%v(i,j,k) = real(v_inc(idx), kind=4)
          if (associated(grid%xa%w)) grid%xa%w(i,j,k) = real(w_inc(idx), kind=4)
          if (associated(grid%xa%t)) grid%xa%t(i,j,k) = real(t_inc(idx), kind=4)
          if (associated(grid%xa%p)) grid%xa%p(i,j,k) = real(p_inc(idx), kind=4)
          if (associated(grid%xa%q)) grid%xa%q(i,j,k) = real(q_inc(idx), kind=4)
          ! Additional moisture/thermodynamic variables (check size - may be dummy arrays)
          if (associated(grid%xa%qt)) then
            if (size(grid%xa%qt, 1) >= nx .and. size(grid%xa%qt, 2) >= ny .and. size(grid%xa%qt, 3) >= nz) then
              grid%xa%qt(i,j,k) = real(qt_inc(idx), kind=4)
            end if
          end if
          if (associated(grid%xa%rh)) then
            if (size(grid%xa%rh, 1) >= nx .and. size(grid%xa%rh, 2) >= ny .and. size(grid%xa%rh, 3) >= nz) then
              grid%xa%rh(i,j,k) = real(rh_inc(idx), kind=4)
            end if
          end if
          ! Derived fields (check size - may also be conditionally allocated)
          if (associated(grid%xa%rho)) then
            if (size(grid%xa%rho, 1) >= nx .and. size(grid%xa%rho, 2) >= ny .and. size(grid%xa%rho, 3) >= nz) then
              grid%xa%rho(i,j,k) = real(rho_inc(idx), kind=4)
            end if
          end if
          if (associated(grid%xa%geoh)) then
            if (size(grid%xa%geoh, 1) >= nx .and. size(grid%xa%geoh, 2) >= ny .and. size(grid%xa%geoh, 3) >= nz) then
              grid%xa%geoh(i,j,k) = real(geoh_inc(idx), kind=4)
            end if
          end if
          if (associated(grid%xa%wh)) then
            if (size(grid%xa%wh, 1) >= nx .and. size(grid%xa%wh, 2) >= ny .and. size(grid%xa%wh, 3) >= nz) then
              grid%xa%wh(i,j,k) = real(wh_inc(idx), kind=4)
            end if
          end if
          ! Hydrometeor fields (check size - conditionally allocated based on cloud_cv_options)
          if (associated(grid%xa%qcw)) then
            if (size(grid%xa%qcw, 1) >= nx .and. size(grid%xa%qcw, 2) >= ny .and. size(grid%xa%qcw, 3) >= nz) then
              grid%xa%qcw(i,j,k) = real(qcw_inc(idx), kind=4)
            end if
          end if
          if (associated(grid%xa%qrn)) then
            if (size(grid%xa%qrn, 1) >= nx .and. size(grid%xa%qrn, 2) >= ny .and. size(grid%xa%qrn, 3) >= nz) then
              grid%xa%qrn(i,j,k) = real(qrn_inc(idx), kind=4)
            end if
          end if
          if (associated(grid%xa%qci)) then
            if (size(grid%xa%qci, 1) >= nx .and. size(grid%xa%qci, 2) >= ny .and. size(grid%xa%qci, 3) >= nz) then
              grid%xa%qci(i,j,k) = real(qci_inc(idx), kind=4)
            end if
          end if
          if (associated(grid%xa%qsn)) then
            if (size(grid%xa%qsn, 1) >= nx .and. size(grid%xa%qsn, 2) >= ny .and. size(grid%xa%qsn, 3) >= nz) then
              grid%xa%qsn(i,j,k) = real(qsn_inc(idx), kind=4)
            end if
          end if
          if (associated(grid%xa%qgr)) then
            if (size(grid%xa%qgr, 1) >= nx .and. size(grid%xa%qgr, 2) >= ny .and. size(grid%xa%qgr, 3) >= nz) then
              grid%xa%qgr(i,j,k) = real(qgr_inc(idx), kind=4)
            end if
          end if
        end do
      end do
    end do
    
    ! Populate 2D fields from C++ arrays (row-major indexing, check size for conditionally allocated fields)
    do j = 1, ny
      do i = 1, nx
        idx = i + (j-1)*nx
        if (associated(grid%xa%psfc)) then
          if (size(grid%xa%psfc, 1) >= nx .and. size(grid%xa%psfc, 2) >= ny) then
            grid%xa%psfc(i,j) = real(psfc_inc(idx), kind=4)
          end if
        end if
        if (associated(grid%xa%mu)) then
          if (size(grid%xa%mu, 1) >= nx .and. size(grid%xa%mu, 2) >= ny) then
            grid%xa%mu(i,j) = real(mu_inc(idx), kind=4)
          end if
        end if
      end do
    end do
    
  end subroutine wrfda_update_analysis_increments

  ! New function to call da_get_innov_vector directly
  integer(c_int) function wrfda_get_innov_vector(it, ob_ptr, iv_ptr, grid_ptr) bind(C, name="wrfda_get_innov_vector")
    implicit none
    integer(c_int), intent(in) :: it
    type(c_ptr), value :: ob_ptr, iv_ptr, grid_ptr
    type(domain), pointer :: grid
    type(y_type), pointer :: ob
    type(iv_type), pointer :: iv
    type(grid_config_rec_type), pointer :: config_flags
    
    ! Convert C pointers to Fortran pointers
    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(ob_ptr, ob)
    call c_f_pointer(iv_ptr, iv)
    
    ! config_flags is required by WRFDA interface - allocate it
    allocate(config_flags)
    
    ! Initialize QC statistics array (imported from da_control module)
    num_qcstat_conv = 0
    
    ! Call the main WRFDA innovation vector computation routine
    call da_get_innov_vector(it, num_qcstat_conv, ob, iv, grid, config_flags)

    ! Clean up allocated config_flags
    deallocate(config_flags)

    wrfda_get_innov_vector = 0
    
  end function wrfda_get_innov_vector


  ! Helper function to construct WRFDA domain structure from flat arrays
  integer(c_int) function wrfda_construct_domain_from_arrays(nx, ny, nz, u, v, t, q, psfc, ph, phb, hf, hgt, p, pb, lats2d, lons2d) bind(C, name="wrfda_construct_domain_from_arrays")
    implicit none
    integer(c_int), intent(in) :: nx, ny, nz
    real(c_double), intent(in) :: u(*), v(*), t(*), q(*), psfc(*)
    real(c_double), intent(in) :: ph(*), phb(*), hf(*), hgt(*), p(*), pb(*)
    real(c_double), intent(in) :: lats2d(*), lons2d(*)
    
    type(domain), pointer :: grid
    integer :: i, j, k, idx, idx_3d
    integer :: staggered_nz  ! Height field has nz+1 vertical levels

    ! Allocate new domain structure
    if (.not.grid_initialized) then
      allocate(persistent_grid)
      grid => persistent_grid
      grid_initialized = .true.
    else
      if (associated(persistent_grid)) then
        grid => persistent_grid
      else
        allocate(persistent_grid)
        grid => persistent_grid
        grid_initialized = .true.
      end if
    end if

    ! Temporary constant for METADA - should be imported from WRFDA constants
    Max_StHeight_Diff = 100.0

    ! Set up basic domain dimensions
    grid%id = 1
    
    ! Initialize domain clock (required for WRFDA time management)
    ! Set domain_clock_created to false initially
    grid%domain_clock_created = .false.
    
    ! Height field is vertically staggered with nz+1 levels
    staggered_nz = nz + 1
    
    ! Allocate PH and PHB fields (vertically staggered with nz+1 levels)
    
    allocate(grid%ph_2(1:nx, 1:ny, 1:staggered_nz))
    allocate(grid%phb(1:nx, 1:ny, 1:staggered_nz))
    

    
    ! Set up grid processor dimensions (required for da_copy_tile_dims)
    grid%xp%kds = 1        ! Start of vertical domain
    grid%xp%kde = nz       ! End of vertical domain (mass variables)
    grid%xp%ids = 1        ! Start of i domain
    grid%xp%ide = nx       ! End of i domain
    grid%xp%jds = 1        ! Start of j domain
    grid%xp%jde = ny       ! End of j domain
    
    ! Set up tile dimensions (single tile for now)
    grid%num_tiles = 1
    allocate(grid%i_start(1:1), grid%i_end(1:1))
    allocate(grid%j_start(1:1), grid%j_end(1:1))
    grid%i_start(1) = 1
    grid%i_end(1) = nx
    grid%j_start(1) = 1
    grid%j_end(1) = ny
    
    ! Set up WRFDA memory dimensions (sm/em arrays)
    grid%sm31 = 1; grid%em31 = nx
    grid%sm32 = 1; grid%em32 = ny
    grid%sm33 = 1; grid%em33 = nz
    grid%sm31x = 1; grid%em31x = nx
    grid%sm32x = 1; grid%em32x = ny
    grid%sm33x = 1; grid%em33x = nz
    
    ! Initialize domain clock for WRFDA time management
    call initialize_domain_clock(grid)
    grid%sm31y = 1; grid%em31y = nx
    grid%sm32y = 1; grid%em32y = ny
    grid%sm33y = 1; grid%em33y = nz
    
    ! Set up WRFDA domain dimensions (sd/ed arrays) - required by da_copy_dims
    ! In WRF, ed arrays are typically the upper bound + 1, so for nx,ny,nz grid:
    grid%sd31 = 1; grid%ed31 = nx + 1
    grid%sd32 = 1; grid%ed32 = ny + 1  
    grid%sd33 = 1; grid%ed33 = nz + 1
    grid%sp31 = 1; grid%ep31 = nx
    grid%sp32 = 1; grid%ep32 = ny
    grid%sp33 = 1; grid%ep33 = nz

    ! Allocate and populate xb (background state) arrays
    allocate(grid%xb%u(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%v(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%t(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%q(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%p(1:nx, 1:ny, 1:nz))
    allocate(grid%xb%psfc(1:nx, 1:ny))

    ! Copy data from flat arrays to WRFDA structure
    ! Note: C++ arrays are in [X,Y,Z] order (row-major), Fortran arrays are [Z,Y,X] (column-major)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx = i + (j-1)*nx + (k-1)*nx*ny
          grid%xb%u(i,j,k) = real(u(idx), kind=4)
          grid%xb%v(i,j,k) = real(v(idx), kind=4)
          grid%xb%t(i,j,k) = real(t(idx), kind=4)
          grid%xb%q(i,j,k) = real(q(idx), kind=4)
          ! Calculate pressure from P and PB: grid%xb%p = pb + p
          grid%xb%p(i,j,k) = real(pb(idx) + p(idx), kind=4)
        end do
      end do
    end do
    
    ! Copy PH and PHB data (vertically staggered with nz+1 levels)
    do k = 1, staggered_nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          grid%ph_2(i,j,k) = real(ph(idx_3d), kind=4)
          grid%phb(i,j,k) = real(phb(idx_3d), kind=4)
        end do
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        ! C++ row-major indexing: [i][j] -> i + j*nx
        idx = i + (j-1)*nx
        grid%xb%psfc(i,j) = real(psfc(idx), kind=4)
      end do
    end do

    ! Set up grid metadata
    allocate(grid%xb%lat(1:nx, 1:ny))
    allocate(grid%xb%lon(1:nx, 1:ny))
    ! Height field (HF) is vertically staggered with nz+1 levels
    allocate(grid%xb%h(1:nx, 1:ny, 1:staggered_nz))
    ! Terrain height fields (HGT) - 2D surface fields
    allocate(grid%ht(1:nx, 1:ny))
    allocate(grid%xb%terr(1:nx, 1:ny))
    ! Surface roughness length - temporary constant value
    allocate(grid%xb%rough(1:nx, 1:ny))

    do j = 1, ny
      do i = 1, nx
        ! C++ row-major indexing: [i][j] -> i + j*nx
        idx = i + (j-1)*nx
        grid%xb%lat(i,j) = real(lats2d(idx), kind=4)
        grid%xb%lon(i,j) = real(lons2d(idx), kind=4)
        ! Assign terrain height from HGT field
        grid%ht(i,j) = real(hgt(idx), kind=4)
        grid%xb%terr(i,j) = real(hgt(idx), kind=4)
        ! Assign temporary constant roughness length
        grid%xb%rough(i,j) = 0.5
        do k = 1, staggered_nz
          ! Use calculated height field instead of levels array
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          ! Note: HF is vertically staggered with nz+1 levels
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          grid%xb%h(i,j,k) = real(hf(idx_3d), kind=4)
        end do
      end do
    end do

    ! Allocate xa fields (analysis increments) and initialize to zero
    ! For adjoint operator, xa starts at zero and gradients are accumulated in jo_grad_x
    allocate(grid%xa%u(1:nx,1:ny,1:nz)); grid%xa%u = 0.0
    allocate(grid%xa%v(1:nx,1:ny,1:nz)); grid%xa%v = 0.0
    allocate(grid%xa%t(1:nx,1:ny,1:nz)); grid%xa%t = 0.0
    allocate(grid%xa%q(1:nx,1:ny,1:nz)); grid%xa%q = 0.0
    allocate(grid%xa%psfc(1:nx,1:ny)); grid%xa%psfc = 0.0

    wrfda_construct_domain_from_arrays = 0
    
  end function wrfda_construct_domain_from_arrays

  ! Proper WRFDA domain initialization using da_transfer_wrftoxb
  ! This function populates the WRF grid structure from arrays and then calls
  ! da_transfer_wrftoxb to properly compute all derived fields
  integer(c_int) function wrfda_init_domain_from_wrf_fields(nx, ny, nz, &
      u, v, w, t, mu, mub, p, pb, ph, phb, &
      xlat, xlong, ht, znu, znw, dn, dnw, rdnw, rdn, &
      p_top, t_init, moist, num_moist, psfc, &
      start_year, start_month, start_day, start_hour) &
      bind(C, name="wrfda_init_domain_from_wrf_fields")
    use da_transfer_model, only: da_transfer_wrftoxb
    use da_define_structures, only: xbx_type
    implicit none
    
    ! Grid dimensions
    integer(c_int), intent(in) :: nx, ny, nz, num_moist
    
    ! WRF state variables (3D fields: nx x ny x nz)
    real(c_double), intent(in) :: u(*), v(*), w(*), t(*)
    real(c_double), intent(in) :: mu(*), mub(*)  ! 2D: nx x ny
    real(c_double), intent(in) :: p(*), pb(*)    ! 3D: nx x ny x nz
    real(c_double), intent(in) :: ph(*), phb(*)  ! 3D staggered: nx x ny x (nz+1)
    
    ! Grid metadata (2D fields: nx x ny)
    real(c_double), intent(in) :: xlat(*), xlong(*), ht(*)
    
    ! Vertical coordinates (1D arrays: size nz or nz+1)
    real(c_double), intent(in) :: znu(*), znw(*), dn(*), dnw(*)
    real(c_double), intent(in) :: rdnw(*), rdn(*)
    
    ! Scalar parameters
    real(c_double), intent(in) :: p_top
    real(c_double), intent(in) :: t_init(*)  ! 3D: nx x ny x nz
    real(c_double), intent(in) :: moist(*)   ! 4D: nx x ny x nz x num_moist
    real(c_double), intent(in) :: psfc(*)    ! 2D: nx x ny
    
    ! Time information
    integer(c_int), intent(in) :: start_year, start_month, start_day, start_hour
    
    ! Local variables
    type(domain), pointer :: grid
    type(xbx_type) :: xbx
    type(grid_config_rec_type) :: config_flags
    integer :: i, j, k, m, idx, idx_3d, staggered_nz
    
    ! Set up grid structure (similar to wrfda_construct_domain_from_arrays but more complete)
    if (.not. grid_initialized) then
      allocate(persistent_grid)
      grid => persistent_grid

      ! CRITICAL: Nullify all pointer components immediately after allocation
      nullify(grid%u_2, grid%v_2, grid%w_2, grid%t_2, grid%p, grid%pb)
      nullify(grid%t_init, grid%mu_2, grid%mub, grid%psfc)
      nullify(grid%ph_2, grid%phb, grid%moist)
      nullify(grid%xlat, grid%xlong, grid%ht)
      nullify(grid%znu, grid%znw, grid%dn, grid%dnw, grid%rdnw, grid%rdn)

      grid_initialized = .true.
    else
      grid => persistent_grid
    end if
    
    staggered_nz = nz + 1
    
    ! Allocate all WRF fields needed by da_transfer_wrftoxb
    ! Note: We only allocate if not already allocated by WRFDA
    ! WRFDA's da_transfer_wrftoxb expects these fields to exist
    
    ! Mass point 3D fields (only allocate if not associated)
    if (.not. associated(grid%u_2)) allocate(grid%u_2(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%v_2)) allocate(grid%v_2(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%w_2)) allocate(grid%w_2(1:nx, 1:ny, 1:staggered_nz))
    
    if (.not. associated(grid%t_2)) allocate(grid%t_2(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%p)) allocate(grid%p(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%pb)) allocate(grid%pb(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%t_init)) allocate(grid%t_init(1:nx, 1:ny, 1:nz))
    
    if (.not. associated(grid%mu_2)) allocate(grid%mu_2(1:nx, 1:ny))
    
    if (.not. associated(grid%mub)) allocate(grid%mub(1:nx, 1:ny))
    
    if (.not. associated(grid%psfc)) allocate(grid%psfc(1:nx, 1:ny))
    
    ! Height fields (staggered, only allocate if not associated)
    if (.not. associated(grid%ph_2)) allocate(grid%ph_2(1:nx, 1:ny, 1:staggered_nz))
    
    if (.not. associated(grid%phb)) allocate(grid%phb(1:nx, 1:ny, 1:staggered_nz))
    
    ! Moisture fields (all species: qvapor, qcloud, qrain, qice, qsnow, qgraupel, etc.)
    if (.not. associated(grid%moist)) allocate(grid%moist(1:nx, 1:ny, 1:nz, 1:num_moist))
    
    ! Grid metadata (only allocate if not associated)
    if (.not. associated(grid%xlat)) allocate(grid%xlat(1:nx, 1:ny))
    
    if (.not. associated(grid%xlong)) allocate(grid%xlong(1:nx, 1:ny))
    
    if (.not. associated(grid%ht)) allocate(grid%ht(1:nx, 1:ny))
    
    ! Vertical coordinates (only allocate if not associated)
    if (.not. associated(grid%znu)) allocate(grid%znu(1:staggered_nz))
    
    if (.not. associated(grid%znw)) allocate(grid%znw(1:staggered_nz))
    
    if (.not. associated(grid%dn)) allocate(grid%dn(1:staggered_nz))
    
    if (.not. associated(grid%dnw)) allocate(grid%dnw(1:staggered_nz))
    
    if (.not. associated(grid%rdnw)) allocate(grid%rdnw(1:staggered_nz))
    
    if (.not. associated(grid%rdn)) allocate(grid%rdn(1:staggered_nz))
    
    ! Copy data from flat C arrays to WRFDA grid structure
    ! 3D fields
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          idx = i + (j-1)*nx + (k-1)*nx*ny
          grid%u_2(i,j,k) = real(u(idx), kind=4)
          grid%v_2(i,j,k) = real(v(idx), kind=4)
          grid%t_2(i,j,k) = real(t(idx), kind=4)
          grid%p(i,j,k) = real(p(idx), kind=4)
          grid%pb(i,j,k) = real(pb(idx), kind=4)
          grid%t_init(i,j,k) = real(t_init(idx), kind=4)
        end do
      end do
    end do
    
    ! Copy moisture fields (all species)
    do m = 1, num_moist
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            ! 4D array indexing: i + (j-1)*nx + (k-1)*nx*ny + (m-1)*nx*ny*nz
            idx = i + (j-1)*nx + (k-1)*nx*ny + (m-1)*nx*ny*nz
            grid%moist(i,j,k,m) = real(moist(idx), kind=4)
          end do
        end do
      end do
    end do
    
    ! W field (staggered)
    do k = 1, staggered_nz
      do j = 1, ny
        do i = 1, nx
          idx_3d = i + (j-1)*nx + (k-1)*nx*ny
          grid%w_2(i,j,k) = real(w(idx_3d), kind=4)
          grid%ph_2(i,j,k) = real(ph(idx_3d), kind=4)
          grid%phb(i,j,k) = real(phb(idx_3d), kind=4)
        end do
      end do
    end do
    
    ! 2D fields
    do j = 1, ny
      do i = 1, nx
        idx = i + (j-1)*nx
        grid%mu_2(i,j) = real(mu(idx), kind=4)
        grid%mub(i,j) = real(mub(idx), kind=4)
        grid%psfc(i,j) = real(psfc(idx), kind=4)
        grid%xlat(i,j) = real(xlat(idx), kind=4)
        grid%xlong(i,j) = real(xlong(idx), kind=4)
        grid%ht(i,j) = real(ht(idx), kind=4)
      end do
    end do
    
    ! 1D vertical coordinates
    do k = 1, nz
      grid%znu(k) = real(znu(k), kind=4)
      grid%znw(k) = real(znw(k), kind=4)
      grid%dn(k) = real(dn(k), kind=4)
      grid%dnw(k) = real(dnw(k), kind=4)
      grid%rdnw(k) = real(rdnw(k), kind=4)
      grid%rdn(k) = real(rdn(k), kind=4)
    end do
    ! Last level for staggered arrays
    grid%znw(staggered_nz) = real(znw(staggered_nz), kind=4)
    grid%dnw(staggered_nz) = real(dnw(staggered_nz), kind=4)
    grid%rdnw(staggered_nz) = real(rdnw(staggered_nz), kind=4)
    
    ! Set scalar parameters
    grid%p_top = real(p_top, kind=4)
    grid%start_year = start_year
    grid%start_month = start_month
    grid%start_day = start_day
    grid%start_hour = start_hour
    
    ! Set up grid dimensions
    grid%xp%ids = 1; grid%xp%ide = nx
    grid%xp%jds = 1; grid%xp%jde = ny
    grid%xp%kds = 1; grid%xp%kde = nz
    
    ! Call da_transfer_wrftoxb to do all the proper initialization
    ! This computes all derived fields, diagnostics, etc.
    call da_transfer_wrftoxb(xbx, grid, config_flags)
    
    wrfda_init_domain_from_wrf_fields = 0  ! Success
    
  end function wrfda_init_domain_from_wrf_fields

  ! Count innovation values from iv_type structure for all observation types
  integer(c_int) function wrfda_count_innovations(iv_ptr, num_innovations) bind(C, name="wrfda_count_innovations")
    implicit none
    type(c_ptr), value :: iv_ptr
    integer(c_int), intent(out) :: num_innovations
    
    type(iv_type), pointer :: iv
    integer :: i, n, count
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(iv_ptr, iv)
    
    if (.not. associated(iv)) then
      num_innovations = 0
      wrfda_count_innovations = 1
      return
    end if
    
    count = 0
    
    ! Count innovations for ALL observation families
    ! Each family has its own observation structure with variables that have .inv field
    
    ! Synop observations
    if (associated(iv%synop) .and. iv%info(synop)%nlocal > 0) then
      do n = 1, iv%info(synop)%nlocal
        if (iv%synop(n)%u%qc == 0 .and. abs(iv%synop(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%synop(n)%v%qc == 0 .and. abs(iv%synop(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%synop(n)%t%qc == 0 .and. abs(iv%synop(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%synop(n)%p%qc == 0 .and. abs(iv%synop(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%synop(n)%q%qc == 0 .and. abs(iv%synop(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Metar observations
    if (associated(iv%metar) .and. iv%info(metar)%nlocal > 0) then
      do n = 1, iv%info(metar)%nlocal
        if (iv%metar(n)%u%qc == 0 .and. abs(iv%metar(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%metar(n)%v%qc == 0 .and. abs(iv%metar(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%metar(n)%t%qc == 0 .and. abs(iv%metar(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%metar(n)%p%qc == 0 .and. abs(iv%metar(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%metar(n)%q%qc == 0 .and. abs(iv%metar(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Ships observations
    if (associated(iv%ships) .and. iv%info(ships)%nlocal > 0) then
      do n = 1, iv%info(ships)%nlocal
        if (iv%ships(n)%u%qc == 0 .and. abs(iv%ships(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%ships(n)%v%qc == 0 .and. abs(iv%ships(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%ships(n)%t%qc == 0 .and. abs(iv%ships(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%ships(n)%p%qc == 0 .and. abs(iv%ships(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%ships(n)%q%qc == 0 .and. abs(iv%ships(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Buoy observations
    if (associated(iv%buoy) .and. iv%info(buoy)%nlocal > 0) then
      do n = 1, iv%info(buoy)%nlocal
        if (iv%buoy(n)%u%qc == 0 .and. abs(iv%buoy(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%buoy(n)%v%qc == 0 .and. abs(iv%buoy(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%buoy(n)%t%qc == 0 .and. abs(iv%buoy(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%buoy(n)%p%qc == 0 .and. abs(iv%buoy(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%buoy(n)%q%qc == 0 .and. abs(iv%buoy(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Sound (radiosonde) observations - multi-level
    if (associated(iv%sound) .and. iv%info(sound)%nlocal > 0) then
      do n = 1, iv%info(sound)%nlocal
        do i = 1, iv%info(sound)%levels(n)
          if (iv%sound(n)%u(i)%qc == 0 .and. abs(iv%sound(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%sound(n)%v(i)%qc == 0 .and. abs(iv%sound(n)%v(i)%inv) > 1.0e-10) count = count + 1
          if (iv%sound(n)%t(i)%qc == 0 .and. abs(iv%sound(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%sound(n)%q(i)%qc == 0 .and. abs(iv%sound(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Airep (aircraft) observations - multi-level
    if (associated(iv%airep) .and. iv%info(airep)%nlocal > 0) then
      do n = 1, iv%info(airep)%nlocal
        do i = 1, iv%info(airep)%levels(n)
          if (iv%airep(n)%u(i)%qc == 0 .and. abs(iv%airep(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%airep(n)%v(i)%qc == 0 .and. abs(iv%airep(n)%v(i)%inv) > 1.0e-10) count = count + 1
          if (iv%airep(n)%t(i)%qc == 0 .and. abs(iv%airep(n)%t(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Pilot (balloon) observations - multi-level
    if (associated(iv%pilot) .and. iv%info(pilot)%nlocal > 0) then
      do n = 1, iv%info(pilot)%nlocal
        do i = 1, iv%info(pilot)%levels(n)
          if (iv%pilot(n)%u(i)%qc == 0 .and. abs(iv%pilot(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%pilot(n)%v(i)%qc == 0 .and. abs(iv%pilot(n)%v(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Sonde surface observations
    if (associated(iv%sonde_sfc) .and. iv%info(sonde_sfc)%nlocal > 0) then
      do n = 1, iv%info(sonde_sfc)%nlocal
        if (iv%sonde_sfc(n)%u%qc == 0 .and. abs(iv%sonde_sfc(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%sonde_sfc(n)%v%qc == 0 .and. abs(iv%sonde_sfc(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%sonde_sfc(n)%t%qc == 0 .and. abs(iv%sonde_sfc(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%sonde_sfc(n)%p%qc == 0 .and. abs(iv%sonde_sfc(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%sonde_sfc(n)%q%qc == 0 .and. abs(iv%sonde_sfc(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Geostationary satellite AMV observations - multi-level
    if (associated(iv%geoamv) .and. iv%info(geoamv)%nlocal > 0) then
      do n = 1, iv%info(geoamv)%nlocal
        do i = 1, iv%info(geoamv)%levels(n)
          if (iv%geoamv(n)%u(i)%qc == 0 .and. abs(iv%geoamv(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%geoamv(n)%v(i)%qc == 0 .and. abs(iv%geoamv(n)%v(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Polar satellite AMV observations - multi-level
    if (associated(iv%polaramv) .and. iv%info(polaramv)%nlocal > 0) then
      do n = 1, iv%info(polaramv)%nlocal
        do i = 1, iv%info(polaramv)%levels(n)
          if (iv%polaramv(n)%u(i)%qc == 0 .and. abs(iv%polaramv(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%polaramv(n)%v(i)%qc == 0 .and. abs(iv%polaramv(n)%v(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! GPS precipitable water observations
    if (associated(iv%gpspw) .and. iv%info(gpspw)%nlocal > 0) then
      do n = 1, iv%info(gpspw)%nlocal
        if (iv%gpspw(n)%tpw%qc == 0 .and. abs(iv%gpspw(n)%tpw%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! GPS refractivity observations - multi-level
    if (associated(iv%gpsref) .and. iv%info(gpsref)%nlocal > 0) then
      do n = 1, iv%info(gpsref)%nlocal
        do i = 1, iv%info(gpsref)%levels(n)
          if (iv%gpsref(n)%ref(i)%qc == 0 .and. abs(iv%gpsref(n)%ref(i)%inv) > 1.0e-10) count = count + 1
          if (iv%gpsref(n)%p(i)%qc == 0 .and. abs(iv%gpsref(n)%p(i)%inv) > 1.0e-10) count = count + 1
          if (iv%gpsref(n)%t(i)%qc == 0 .and. abs(iv%gpsref(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%gpsref(n)%q(i)%qc == 0 .and. abs(iv%gpsref(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! GPS excess phase observations - multi-level
    if (associated(iv%gpseph) .and. iv%info(gpseph)%nlocal > 0) then
      do n = 1, iv%info(gpseph)%nlocal
        do i = 1, iv%info(gpseph)%levels(n)
          if (iv%gpseph(n)%eph(i)%qc == 0 .and. abs(iv%gpseph(n)%eph(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! QuikSCAT wind observations
    if (associated(iv%qscat) .and. iv%info(qscat)%nlocal > 0) then
      do n = 1, iv%info(qscat)%nlocal
        if (iv%qscat(n)%u%qc == 0 .and. abs(iv%qscat(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%qscat(n)%v%qc == 0 .and. abs(iv%qscat(n)%v%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Wind profiler observations - multi-level (uses pilot_type)
    if (associated(iv%profiler) .and. iv%info(profiler)%nlocal > 0) then
      do n = 1, iv%info(profiler)%nlocal
        do i = 1, iv%info(profiler)%levels(n)
          if (iv%profiler(n)%u(i)%qc == 0 .and. abs(iv%profiler(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%profiler(n)%v(i)%qc == 0 .and. abs(iv%profiler(n)%v(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! SSM/I rain rate observations
    if (associated(iv%ssmi_rv) .and. iv%info(ssmi_rv)%nlocal > 0) then
      do n = 1, iv%info(ssmi_rv)%nlocal
        if (iv%ssmi_rv(n)%Speed%qc == 0 .and. abs(iv%ssmi_rv(n)%Speed%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_rv(n)%tpw%qc == 0 .and. abs(iv%ssmi_rv(n)%tpw%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! SSM/I brightness temperature observations - single level, multiple channels
    if (associated(iv%ssmi_tb) .and. iv%info(ssmi_tb)%nlocal > 0) then
      do n = 1, iv%info(ssmi_tb)%nlocal
        if (iv%ssmi_tb(n)%tb19h%qc == 0 .and. abs(iv%ssmi_tb(n)%tb19h%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb19v%qc == 0 .and. abs(iv%ssmi_tb(n)%tb19v%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb22v%qc == 0 .and. abs(iv%ssmi_tb(n)%tb22v%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb37h%qc == 0 .and. abs(iv%ssmi_tb(n)%tb37h%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb37v%qc == 0 .and. abs(iv%ssmi_tb(n)%tb37v%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb85h%qc == 0 .and. abs(iv%ssmi_tb(n)%tb85h%inv) > 1.0e-10) count = count + 1
        if (iv%ssmi_tb(n)%tb85v%qc == 0 .and. abs(iv%ssmi_tb(n)%tb85v%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! SSM/T1 observations - multi-level
    if (associated(iv%ssmt1) .and. iv%info(ssmt1)%nlocal > 0) then
      do n = 1, iv%info(ssmt1)%nlocal
        do i = 1, iv%info(ssmt1)%levels(n)
          if (iv%ssmt1(n)%t(i)%qc == 0 .and. abs(iv%ssmt1(n)%t(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! SSM/T2 observations - multi-level
    if (associated(iv%ssmt2) .and. iv%info(ssmt2)%nlocal > 0) then
      do n = 1, iv%info(ssmt2)%nlocal
        do i = 1, iv%info(ssmt2)%levels(n)
          if (iv%ssmt2(n)%rh(i)%qc == 0 .and. abs(iv%ssmt2(n)%rh(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! SATEM observations - multi-level
    if (associated(iv%satem) .and. iv%info(satem)%nlocal > 0) then
      do n = 1, iv%info(satem)%nlocal
        do i = 1, iv%info(satem)%levels(n)
          if (iv%satem(n)%thickness(i)%qc == 0 .and. abs(iv%satem(n)%thickness(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Pseudo observations - single level
    if (associated(iv%pseudo) .and. iv%info(pseudo)%nlocal > 0) then
      do n = 1, iv%info(pseudo)%nlocal
        if (iv%pseudo(n)%u%qc == 0 .and. abs(iv%pseudo(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%pseudo(n)%v%qc == 0 .and. abs(iv%pseudo(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%pseudo(n)%t%qc == 0 .and. abs(iv%pseudo(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%pseudo(n)%p%qc == 0 .and. abs(iv%pseudo(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%pseudo(n)%q%qc == 0 .and. abs(iv%pseudo(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Bogus (tropical cyclone) observations - multi-level
    if (associated(iv%bogus) .and. iv%info(bogus)%nlocal > 0) then
      do n = 1, iv%info(bogus)%nlocal
        do i = 1, iv%info(bogus)%levels(n)
          if (iv%bogus(n)%u(i)%qc == 0 .and. abs(iv%bogus(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%bogus(n)%v(i)%qc == 0 .and. abs(iv%bogus(n)%v(i)%inv) > 1.0e-10) count = count + 1
          if (iv%bogus(n)%t(i)%qc == 0 .and. abs(iv%bogus(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%bogus(n)%q(i)%qc == 0 .and. abs(iv%bogus(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
        ! Surface level
        if (iv%bogus(n)%slp%qc == 0 .and. abs(iv%bogus(n)%slp%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! AIRS retrievals - multi-level
    if (associated(iv%airsr) .and. iv%info(airsr)%nlocal > 0) then
      do n = 1, iv%info(airsr)%nlocal
        do i = 1, iv%info(airsr)%levels(n)
          if (iv%airsr(n)%t(i)%qc == 0 .and. abs(iv%airsr(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%airsr(n)%q(i)%qc == 0 .and. abs(iv%airsr(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! MTGIRS observations - multi-level
    if (associated(iv%mtgirs) .and. iv%info(mtgirs)%nlocal > 0) then
      do n = 1, iv%info(mtgirs)%nlocal
        do i = 1, iv%info(mtgirs)%levels(n)
          if (iv%mtgirs(n)%u(i)%qc == 0 .and. abs(iv%mtgirs(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%mtgirs(n)%v(i)%qc == 0 .and. abs(iv%mtgirs(n)%v(i)%inv) > 1.0e-10) count = count + 1
          if (iv%mtgirs(n)%t(i)%qc == 0 .and. abs(iv%mtgirs(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%mtgirs(n)%q(i)%qc == 0 .and. abs(iv%mtgirs(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! TAMDAR aircraft observations - multi-level
    if (associated(iv%tamdar) .and. iv%info(tamdar)%nlocal > 0) then
      do n = 1, iv%info(tamdar)%nlocal
        do i = 1, iv%info(tamdar)%levels(n)
          if (iv%tamdar(n)%u(i)%qc == 0 .and. abs(iv%tamdar(n)%u(i)%inv) > 1.0e-10) count = count + 1
          if (iv%tamdar(n)%v(i)%qc == 0 .and. abs(iv%tamdar(n)%v(i)%inv) > 1.0e-10) count = count + 1
          if (iv%tamdar(n)%t(i)%qc == 0 .and. abs(iv%tamdar(n)%t(i)%inv) > 1.0e-10) count = count + 1
          if (iv%tamdar(n)%q(i)%qc == 0 .and. abs(iv%tamdar(n)%q(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! TAMDAR surface observations (uses synop_type)
    if (associated(iv%tamdar_sfc) .and. iv%info(tamdar_sfc)%nlocal > 0) then
      do n = 1, iv%info(tamdar_sfc)%nlocal
        if (iv%tamdar_sfc(n)%u%qc == 0 .and. abs(iv%tamdar_sfc(n)%u%inv) > 1.0e-10) count = count + 1
        if (iv%tamdar_sfc(n)%v%qc == 0 .and. abs(iv%tamdar_sfc(n)%v%inv) > 1.0e-10) count = count + 1
        if (iv%tamdar_sfc(n)%t%qc == 0 .and. abs(iv%tamdar_sfc(n)%t%inv) > 1.0e-10) count = count + 1
        if (iv%tamdar_sfc(n)%p%qc == 0 .and. abs(iv%tamdar_sfc(n)%p%inv) > 1.0e-10) count = count + 1
        if (iv%tamdar_sfc(n)%q%qc == 0 .and. abs(iv%tamdar_sfc(n)%q%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Rain gauge observations
    if (associated(iv%rain) .and. iv%info(rain)%nlocal > 0) then
      do n = 1, iv%info(rain)%nlocal
        if (iv%rain(n)%rain%qc == 0 .and. abs(iv%rain(n)%rain%inv) > 1.0e-10) count = count + 1
      end do
    end if
    
    ! Radar observations - multi-level
    if (associated(iv%radar) .and. iv%info(radar)%nlocal > 0) then
      do n = 1, iv%info(radar)%nlocal
        do i = 1, iv%info(radar)%levels(n)
          if (iv%radar(n)%rv(i)%qc == 0 .and. abs(iv%radar(n)%rv(i)%inv) > 1.0e-10) count = count + 1
          if (iv%radar(n)%rf(i)%qc == 0 .and. abs(iv%radar(n)%rf(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    ! Lightning observations
    if (associated(iv%lightning) .and. iv%info(lightning)%nlocal > 0) then
      do n = 1, iv%info(lightning)%nlocal
        do i = 1, iv%info(lightning)%levels(n)
          if (iv%lightning(n)%w(i)%qc == 0 .and. abs(iv%lightning(n)%w(i)%inv) > 1.0e-10) count = count + 1
          if (iv%lightning(n)%div(i)%qc == 0 .and. abs(iv%lightning(n)%div(i)%inv) > 1.0e-10) count = count + 1
          if (iv%lightning(n)%qv(i)%qc == 0 .and. abs(iv%lightning(n)%qv(i)%inv) > 1.0e-10) count = count + 1
        end do
      end do
    end if
    
    num_innovations = count
    wrfda_count_innovations = 0
    
  end function wrfda_count_innovations

  ! Extract innovation values from iv_type structure
  integer(c_int) function wrfda_extract_innovations(iv_ptr, innovations, num_innovations) bind(C, name="wrfda_extract_innovations")
    implicit none
    type(c_ptr), value :: iv_ptr
    real(c_double), intent(out) :: innovations(*)
    integer(c_int), intent(out) :: num_innovations
    
    type(iv_type), pointer :: iv
    integer :: i, k, count
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(iv_ptr, iv)
    
    if (.not. associated(iv)) then
      num_innovations = 0
      wrfda_extract_innovations = 1
      return
    end if
    
    ! Extract innovations for all families
    count = 0
    
    ! SYNOP family
    if (associated(iv%synop) .and. iv%info(synop)%nlocal > 0) then
      do i = 1, iv%info(synop)%nlocal
        if (i > size(iv%synop)) exit
        
        ! Extract U component innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%u%qc == 0 .and. abs(iv%synop(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%synop(i)%u%inv
        end if
        
        ! Extract V component innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%v%qc == 0 .and. abs(iv%synop(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%synop(i)%v%inv
        end if
        
        ! Extract T component innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%t%qc == 0 .and. abs(iv%synop(i)%t%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%synop(i)%t%inv
        end if
        
        ! Extract P component innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%p%qc == 0 .and. abs(iv%synop(i)%p%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%synop(i)%p%inv
        end if
        
        ! Extract Q component innovation if QC is good and innovation is non-zero
        if (iv%synop(i)%q%qc == 0 .and. abs(iv%synop(i)%q%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%synop(i)%q%inv
        end if
      end do
    end if
    
    ! METAR family
    if (associated(iv%metar) .and. iv%info(metar)%nlocal > 0) then
      do i = 1, iv%info(metar)%nlocal
        if (i > size(iv%metar)) exit
        
        if (iv%metar(i)%u%qc == 0 .and. abs(iv%metar(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%metar(i)%u%inv
        end if
        if (iv%metar(i)%v%qc == 0 .and. abs(iv%metar(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%metar(i)%v%inv
        end if
        if (iv%metar(i)%t%qc == 0 .and. abs(iv%metar(i)%t%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%metar(i)%t%inv
        end if
        if (iv%metar(i)%p%qc == 0 .and. abs(iv%metar(i)%p%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%metar(i)%p%inv
        end if
        if (iv%metar(i)%q%qc == 0 .and. abs(iv%metar(i)%q%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%metar(i)%q%inv
        end if
      end do
    end if
    
    ! SHIPS family
    if (associated(iv%ships) .and. iv%info(ships)%nlocal > 0) then
      do i = 1, iv%info(ships)%nlocal
        if (i > size(iv%ships)) exit
        
        if (iv%ships(i)%u%qc == 0 .and. abs(iv%ships(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ships(i)%u%inv
        end if
        if (iv%ships(i)%v%qc == 0 .and. abs(iv%ships(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ships(i)%v%inv
        end if
        if (iv%ships(i)%t%qc == 0 .and. abs(iv%ships(i)%t%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ships(i)%t%inv
        end if
        if (iv%ships(i)%p%qc == 0 .and. abs(iv%ships(i)%p%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ships(i)%p%inv
        end if
        if (iv%ships(i)%q%qc == 0 .and. abs(iv%ships(i)%q%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ships(i)%q%inv
        end if
      end do
    end if
    
    ! SOUND family - if available and allocated
    ! Note: sound observations have multiple levels per station
    if (associated(iv%sound) .and. iv%info(sound)%nlocal > 0) then
      do i = 1, iv%info(sound)%nlocal
        if (i > size(iv%sound)) exit
        
        ! Loop through levels for this station
        do k = 1, iv%info(sound)%levels(i)
          ! Extract U component if QC is good and innovation is non-zero
          if (iv%sound(i)%u(k)%qc == 0 .and. abs(iv%sound(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%sound(i)%u(k)%inv
          end if
          
          ! Extract V component if QC is good and innovation is non-zero
          if (iv%sound(i)%v(k)%qc == 0 .and. abs(iv%sound(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%sound(i)%v(k)%inv
          end if
          
          ! Extract T component if QC is good and innovation is non-zero
          if (iv%sound(i)%t(k)%qc == 0 .and. abs(iv%sound(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%sound(i)%t(k)%inv
          end if
          
          ! Extract Q component if QC is good and innovation is non-zero
          if (iv%sound(i)%q(k)%qc == 0 .and. abs(iv%sound(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%sound(i)%q(k)%inv
          end if
        end do
      end do
    end if
    
    ! BUOY family
    if (associated(iv%buoy) .and. iv%info(buoy)%nlocal > 0) then
      do i = 1, iv%info(buoy)%nlocal
        if (i > size(iv%buoy)) exit
        
        if (iv%buoy(i)%u%qc == 0 .and. abs(iv%buoy(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%buoy(i)%u%inv
        end if
        if (iv%buoy(i)%v%qc == 0 .and. abs(iv%buoy(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%buoy(i)%v%inv
        end if
        if (iv%buoy(i)%t%qc == 0 .and. abs(iv%buoy(i)%t%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%buoy(i)%t%inv
        end if
        if (iv%buoy(i)%p%qc == 0 .and. abs(iv%buoy(i)%p%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%buoy(i)%p%inv
        end if
        if (iv%buoy(i)%q%qc == 0 .and. abs(iv%buoy(i)%q%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%buoy(i)%q%inv
        end if
      end do
    end if
    
    ! AIREP family - multi-level
    if (associated(iv%airep) .and. iv%info(airep)%nlocal > 0) then
      do i = 1, iv%info(airep)%nlocal
        if (i > size(iv%airep)) exit
        
        do k = 1, iv%info(airep)%levels(i)
          if (iv%airep(i)%u(k)%qc == 0 .and. abs(iv%airep(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%airep(i)%u(k)%inv
          end if
          if (iv%airep(i)%v(k)%qc == 0 .and. abs(iv%airep(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%airep(i)%v(k)%inv
          end if
          if (iv%airep(i)%t(k)%qc == 0 .and. abs(iv%airep(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%airep(i)%t(k)%inv
          end if
        end do
      end do
    end if
    
    ! PILOT family - multi-level
    if (associated(iv%pilot) .and. iv%info(pilot)%nlocal > 0) then
      do i = 1, iv%info(pilot)%nlocal
        if (i > size(iv%pilot)) exit
        
        do k = 1, iv%info(pilot)%levels(i)
          if (iv%pilot(i)%u(k)%qc == 0 .and. abs(iv%pilot(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%pilot(i)%u(k)%inv
          end if
          if (iv%pilot(i)%v(k)%qc == 0 .and. abs(iv%pilot(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%pilot(i)%v(k)%inv
          end if
        end do
      end do
    end if
    
    ! SONDE_SFC family
    if (associated(iv%sonde_sfc) .and. iv%info(sonde_sfc)%nlocal > 0) then
      do i = 1, iv%info(sonde_sfc)%nlocal
        if (i > size(iv%sonde_sfc)) exit
        
        if (iv%sonde_sfc(i)%u%qc == 0 .and. abs(iv%sonde_sfc(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%sonde_sfc(i)%u%inv
        end if
        if (iv%sonde_sfc(i)%v%qc == 0 .and. abs(iv%sonde_sfc(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%sonde_sfc(i)%v%inv
        end if
        if (iv%sonde_sfc(i)%t%qc == 0 .and. abs(iv%sonde_sfc(i)%t%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%sonde_sfc(i)%t%inv
        end if
        if (iv%sonde_sfc(i)%p%qc == 0 .and. abs(iv%sonde_sfc(i)%p%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%sonde_sfc(i)%p%inv
        end if
        if (iv%sonde_sfc(i)%q%qc == 0 .and. abs(iv%sonde_sfc(i)%q%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%sonde_sfc(i)%q%inv
        end if
      end do
    end if
    
    ! GEOAMV family - multi-level
    if (associated(iv%geoamv) .and. iv%info(geoamv)%nlocal > 0) then
      do i = 1, iv%info(geoamv)%nlocal
        if (i > size(iv%geoamv)) exit
        
        do k = 1, iv%info(geoamv)%levels(i)
          if (iv%geoamv(i)%u(k)%qc == 0 .and. abs(iv%geoamv(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%geoamv(i)%u(k)%inv
          end if
          if (iv%geoamv(i)%v(k)%qc == 0 .and. abs(iv%geoamv(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%geoamv(i)%v(k)%inv
          end if
        end do
      end do
    end if
    
    ! POLARAMV family - multi-level
    if (associated(iv%polaramv) .and. iv%info(polaramv)%nlocal > 0) then
      do i = 1, iv%info(polaramv)%nlocal
        if (i > size(iv%polaramv)) exit
        
        do k = 1, iv%info(polaramv)%levels(i)
          if (iv%polaramv(i)%u(k)%qc == 0 .and. abs(iv%polaramv(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%polaramv(i)%u(k)%inv
          end if
          if (iv%polaramv(i)%v(k)%qc == 0 .and. abs(iv%polaramv(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%polaramv(i)%v(k)%inv
          end if
        end do
      end do
    end if
    
    ! GPSPW family
    if (associated(iv%gpspw) .and. iv%info(gpspw)%nlocal > 0) then
      do i = 1, iv%info(gpspw)%nlocal
        if (i > size(iv%gpspw)) exit
        
        if (iv%gpspw(i)%tpw%qc == 0 .and. abs(iv%gpspw(i)%tpw%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%gpspw(i)%tpw%inv
        end if
      end do
    end if
    
    ! GPSREF family - multi-level
    if (associated(iv%gpsref) .and. iv%info(gpsref)%nlocal > 0) then
      do i = 1, iv%info(gpsref)%nlocal
        if (i > size(iv%gpsref)) exit
        
        do k = 1, iv%info(gpsref)%levels(i)
          if (iv%gpsref(i)%ref(k)%qc == 0 .and. abs(iv%gpsref(i)%ref(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%gpsref(i)%ref(k)%inv
          end if
          if (iv%gpsref(i)%p(k)%qc == 0 .and. abs(iv%gpsref(i)%p(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%gpsref(i)%p(k)%inv
          end if
          if (iv%gpsref(i)%t(k)%qc == 0 .and. abs(iv%gpsref(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%gpsref(i)%t(k)%inv
          end if
          if (iv%gpsref(i)%q(k)%qc == 0 .and. abs(iv%gpsref(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%gpsref(i)%q(k)%inv
          end if
        end do
      end do
    end if
    
    ! GPSEPH family - multi-level
    if (associated(iv%gpseph) .and. iv%info(gpseph)%nlocal > 0) then
      do i = 1, iv%info(gpseph)%nlocal
        if (i > size(iv%gpseph)) exit
        
        do k = 1, iv%info(gpseph)%levels(i)
          if (iv%gpseph(i)%eph(k)%qc == 0 .and. abs(iv%gpseph(i)%eph(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%gpseph(i)%eph(k)%inv
          end if
        end do
      end do
    end if
    
    ! QSCAT family
    if (associated(iv%qscat) .and. iv%info(qscat)%nlocal > 0) then
      do i = 1, iv%info(qscat)%nlocal
        if (i > size(iv%qscat)) exit
        
        if (iv%qscat(i)%u%qc == 0 .and. abs(iv%qscat(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%qscat(i)%u%inv
        end if
        if (iv%qscat(i)%v%qc == 0 .and. abs(iv%qscat(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%qscat(i)%v%inv
        end if
      end do
    end if
    
    ! PROFILER family - multi-level
    if (associated(iv%profiler) .and. iv%info(profiler)%nlocal > 0) then
      do i = 1, iv%info(profiler)%nlocal
        if (i > size(iv%profiler)) exit
        
        do k = 1, iv%info(profiler)%levels(i)
          if (iv%profiler(i)%u(k)%qc == 0 .and. abs(iv%profiler(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%profiler(i)%u(k)%inv
          end if
          if (iv%profiler(i)%v(k)%qc == 0 .and. abs(iv%profiler(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%profiler(i)%v(k)%inv
          end if
        end do
      end do
    end if
    
    ! SSMI_RV family
    if (associated(iv%ssmi_rv) .and. iv%info(ssmi_rv)%nlocal > 0) then
      do i = 1, iv%info(ssmi_rv)%nlocal
        if (i > size(iv%ssmi_rv)) exit
        
        if (iv%ssmi_rv(i)%Speed%qc == 0 .and. abs(iv%ssmi_rv(i)%Speed%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ssmi_rv(i)%Speed%inv
        end if
        if (iv%ssmi_rv(i)%tpw%qc == 0 .and. abs(iv%ssmi_rv(i)%tpw%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ssmi_rv(i)%tpw%inv
        end if
      end do
    end if
    
    ! SSMI_TB family
    if (associated(iv%ssmi_tb) .and. iv%info(ssmi_tb)%nlocal > 0) then
      do i = 1, iv%info(ssmi_tb)%nlocal
        if (i > size(iv%ssmi_tb)) exit
        
        if (iv%ssmi_tb(i)%tb19h%qc == 0 .and. abs(iv%ssmi_tb(i)%tb19h%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ssmi_tb(i)%tb19h%inv
        end if
        if (iv%ssmi_tb(i)%tb19v%qc == 0 .and. abs(iv%ssmi_tb(i)%tb19v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ssmi_tb(i)%tb19v%inv
        end if
        if (iv%ssmi_tb(i)%tb22v%qc == 0 .and. abs(iv%ssmi_tb(i)%tb22v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ssmi_tb(i)%tb22v%inv
        end if
        if (iv%ssmi_tb(i)%tb37h%qc == 0 .and. abs(iv%ssmi_tb(i)%tb37h%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ssmi_tb(i)%tb37h%inv
        end if
        if (iv%ssmi_tb(i)%tb37v%qc == 0 .and. abs(iv%ssmi_tb(i)%tb37v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ssmi_tb(i)%tb37v%inv
        end if
        if (iv%ssmi_tb(i)%tb85h%qc == 0 .and. abs(iv%ssmi_tb(i)%tb85h%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ssmi_tb(i)%tb85h%inv
        end if
        if (iv%ssmi_tb(i)%tb85v%qc == 0 .and. abs(iv%ssmi_tb(i)%tb85v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%ssmi_tb(i)%tb85v%inv
        end if
      end do
    end if
    
    ! SSMT1 family - multi-level
    if (associated(iv%ssmt1) .and. iv%info(ssmt1)%nlocal > 0) then
      do i = 1, iv%info(ssmt1)%nlocal
        if (i > size(iv%ssmt1)) exit
        
        do k = 1, iv%info(ssmt1)%levels(i)
          if (iv%ssmt1(i)%t(k)%qc == 0 .and. abs(iv%ssmt1(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%ssmt1(i)%t(k)%inv
          end if
        end do
      end do
    end if
    
    ! SSMT2 family - multi-level
    if (associated(iv%ssmt2) .and. iv%info(ssmt2)%nlocal > 0) then
      do i = 1, iv%info(ssmt2)%nlocal
        if (i > size(iv%ssmt2)) exit
        
        do k = 1, iv%info(ssmt2)%levels(i)
          if (iv%ssmt2(i)%rh(k)%qc == 0 .and. abs(iv%ssmt2(i)%rh(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%ssmt2(i)%rh(k)%inv
          end if
        end do
      end do
    end if
    
    ! SATEM family - multi-level
    if (associated(iv%satem) .and. iv%info(satem)%nlocal > 0) then
      do i = 1, iv%info(satem)%nlocal
        if (i > size(iv%satem)) exit
        
        do k = 1, iv%info(satem)%levels(i)
          if (iv%satem(i)%thickness(k)%qc == 0 .and. abs(iv%satem(i)%thickness(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%satem(i)%thickness(k)%inv
          end if
        end do
      end do
    end if
    
    ! PSEUDO family
    if (associated(iv%pseudo) .and. iv%info(pseudo)%nlocal > 0) then
      do i = 1, iv%info(pseudo)%nlocal
        if (i > size(iv%pseudo)) exit
        
        if (iv%pseudo(i)%u%qc == 0 .and. abs(iv%pseudo(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%pseudo(i)%u%inv
        end if
        if (iv%pseudo(i)%v%qc == 0 .and. abs(iv%pseudo(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%pseudo(i)%v%inv
        end if
        if (iv%pseudo(i)%t%qc == 0 .and. abs(iv%pseudo(i)%t%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%pseudo(i)%t%inv
        end if
        if (iv%pseudo(i)%p%qc == 0 .and. abs(iv%pseudo(i)%p%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%pseudo(i)%p%inv
        end if
        if (iv%pseudo(i)%q%qc == 0 .and. abs(iv%pseudo(i)%q%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%pseudo(i)%q%inv
        end if
      end do
    end if
    
    ! BOGUS family - multi-level
    if (associated(iv%bogus) .and. iv%info(bogus)%nlocal > 0) then
      do i = 1, iv%info(bogus)%nlocal
        if (i > size(iv%bogus)) exit
        
        do k = 1, iv%info(bogus)%levels(i)
          if (iv%bogus(i)%u(k)%qc == 0 .and. abs(iv%bogus(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%bogus(i)%u(k)%inv
          end if
          if (iv%bogus(i)%v(k)%qc == 0 .and. abs(iv%bogus(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%bogus(i)%v(k)%inv
          end if
          if (iv%bogus(i)%t(k)%qc == 0 .and. abs(iv%bogus(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%bogus(i)%t(k)%inv
          end if
          if (iv%bogus(i)%q(k)%qc == 0 .and. abs(iv%bogus(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%bogus(i)%q(k)%inv
          end if
        end do
        if (iv%bogus(i)%slp%qc == 0 .and. abs(iv%bogus(i)%slp%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%bogus(i)%slp%inv
        end if
      end do
    end if
    
    ! AIRSR family - multi-level
    if (associated(iv%airsr) .and. iv%info(airsr)%nlocal > 0) then
      do i = 1, iv%info(airsr)%nlocal
        if (i > size(iv%airsr)) exit
        
        do k = 1, iv%info(airsr)%levels(i)
          if (iv%airsr(i)%t(k)%qc == 0 .and. abs(iv%airsr(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%airsr(i)%t(k)%inv
          end if
          if (iv%airsr(i)%q(k)%qc == 0 .and. abs(iv%airsr(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%airsr(i)%q(k)%inv
          end if
        end do
      end do
    end if
    
    ! MTGIRS family - multi-level
    if (associated(iv%mtgirs) .and. iv%info(mtgirs)%nlocal > 0) then
      do i = 1, iv%info(mtgirs)%nlocal
        if (i > size(iv%mtgirs)) exit
        
        do k = 1, iv%info(mtgirs)%levels(i)
          if (iv%mtgirs(i)%u(k)%qc == 0 .and. abs(iv%mtgirs(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%mtgirs(i)%u(k)%inv
          end if
          if (iv%mtgirs(i)%v(k)%qc == 0 .and. abs(iv%mtgirs(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%mtgirs(i)%v(k)%inv
          end if
          if (iv%mtgirs(i)%t(k)%qc == 0 .and. abs(iv%mtgirs(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%mtgirs(i)%t(k)%inv
          end if
          if (iv%mtgirs(i)%q(k)%qc == 0 .and. abs(iv%mtgirs(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%mtgirs(i)%q(k)%inv
          end if
        end do
      end do
    end if
    
    ! TAMDAR family - multi-level
    if (associated(iv%tamdar) .and. iv%info(tamdar)%nlocal > 0) then
      do i = 1, iv%info(tamdar)%nlocal
        if (i > size(iv%tamdar)) exit
        
        do k = 1, iv%info(tamdar)%levels(i)
          if (iv%tamdar(i)%u(k)%qc == 0 .and. abs(iv%tamdar(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%tamdar(i)%u(k)%inv
          end if
          if (iv%tamdar(i)%v(k)%qc == 0 .and. abs(iv%tamdar(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%tamdar(i)%v(k)%inv
          end if
          if (iv%tamdar(i)%t(k)%qc == 0 .and. abs(iv%tamdar(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%tamdar(i)%t(k)%inv
          end if
          if (iv%tamdar(i)%q(k)%qc == 0 .and. abs(iv%tamdar(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%tamdar(i)%q(k)%inv
          end if
        end do
      end do
    end if
    
    ! TAMDAR_SFC family
    if (associated(iv%tamdar_sfc) .and. iv%info(tamdar_sfc)%nlocal > 0) then
      do i = 1, iv%info(tamdar_sfc)%nlocal
        if (i > size(iv%tamdar_sfc)) exit
        
        if (iv%tamdar_sfc(i)%u%qc == 0 .and. abs(iv%tamdar_sfc(i)%u%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%tamdar_sfc(i)%u%inv
        end if
        if (iv%tamdar_sfc(i)%v%qc == 0 .and. abs(iv%tamdar_sfc(i)%v%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%tamdar_sfc(i)%v%inv
        end if
        if (iv%tamdar_sfc(i)%t%qc == 0 .and. abs(iv%tamdar_sfc(i)%t%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%tamdar_sfc(i)%t%inv
        end if
        if (iv%tamdar_sfc(i)%p%qc == 0 .and. abs(iv%tamdar_sfc(i)%p%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%tamdar_sfc(i)%p%inv
        end if
        if (iv%tamdar_sfc(i)%q%qc == 0 .and. abs(iv%tamdar_sfc(i)%q%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%tamdar_sfc(i)%q%inv
        end if
      end do
    end if
    
    ! RAIN family
    if (associated(iv%rain) .and. iv%info(rain)%nlocal > 0) then
      do i = 1, iv%info(rain)%nlocal
        if (i > size(iv%rain)) exit
        
        if (iv%rain(i)%rain%qc == 0 .and. abs(iv%rain(i)%rain%inv) > 1.0e-10) then
          count = count + 1
          innovations(count) = iv%rain(i)%rain%inv
        end if
      end do
    end if
    
    ! RADAR family - multi-level
    if (associated(iv%radar) .and. iv%info(radar)%nlocal > 0) then
      do i = 1, iv%info(radar)%nlocal
        if (i > size(iv%radar)) exit
        
        do k = 1, iv%info(radar)%levels(i)
          if (iv%radar(i)%rv(k)%qc == 0 .and. abs(iv%radar(i)%rv(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%radar(i)%rv(k)%inv
          end if
          if (iv%radar(i)%rf(k)%qc == 0 .and. abs(iv%radar(i)%rf(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%radar(i)%rf(k)%inv
          end if
        end do
      end do
    end if
    
    ! LIGHTNING family - multi-level
    if (associated(iv%lightning) .and. iv%info(lightning)%nlocal > 0) then
      do i = 1, iv%info(lightning)%nlocal
        if (i > size(iv%lightning)) exit
        
        do k = 1, iv%info(lightning)%levels(i)
          if (iv%lightning(i)%w(k)%qc == 0 .and. abs(iv%lightning(i)%w(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%lightning(i)%w(k)%inv
          end if
          if (iv%lightning(i)%div(k)%qc == 0 .and. abs(iv%lightning(i)%div(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%lightning(i)%div(k)%inv
          end if
          if (iv%lightning(i)%qv(k)%qc == 0 .and. abs(iv%lightning(i)%qv(k)%inv) > 1.0e-10) then
            count = count + 1
            innovations(count) = iv%lightning(i)%qv(k)%inv
          end if
        end do
      end do
    end if
    
    num_innovations = count
    wrfda_extract_innovations = 0
    
  end function wrfda_extract_innovations

  ! Extract observation values from y_type structure for all families
  integer(c_int) function wrfda_extract_observations(iv_ptr, ob_ptr, observations, num_observations) bind(C, name="wrfda_extract_observations")
    implicit none
    type(c_ptr), value :: iv_ptr
    type(c_ptr), value :: ob_ptr
    real(c_double), intent(out) :: observations(*)
    integer(c_int), intent(out) :: num_observations
    
    type(y_type), pointer :: y
    type(iv_type), pointer :: iv
    integer :: i, k, count
    
    ! Convert C pointers to Fortran pointers
    call c_f_pointer(iv_ptr, iv)
    call c_f_pointer(ob_ptr, y)
    
    if (.not. associated(iv)) then
      num_observations = 0
      wrfda_extract_observations = 1
      return
    end if
    
    if (.not. associated(y)) then
      num_observations = 0
      wrfda_extract_observations = 1
      return
    end if
    
    ! Extract observations for all families
    count = 0
    
    ! SYNOP family
    if (associated(y%synop) .and. associated(iv%synop) .and. iv%info(synop)%nlocal > 0) then
      do i = 1, min(size(y%synop), iv%info(synop)%nlocal)
        ! Extract U component if QC is good and innovation is non-zero
        if (iv%synop(i)%u%qc == 0 .and. abs(iv%synop(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%synop(i)%u
        end if
        
        ! Extract V component if QC is good and innovation is non-zero
        if (iv%synop(i)%v%qc == 0 .and. abs(iv%synop(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%synop(i)%v
        end if
        
        ! Extract T component if QC is good and innovation is non-zero
        if (iv%synop(i)%t%qc == 0 .and. abs(iv%synop(i)%t%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%synop(i)%t
        end if
        
        ! Extract P component if QC is good and innovation is non-zero
        if (iv%synop(i)%p%qc == 0 .and. abs(iv%synop(i)%p%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%synop(i)%p
        end if
        
        ! Extract Q component if QC is good and innovation is non-zero
        if (iv%synop(i)%q%qc == 0 .and. abs(iv%synop(i)%q%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%synop(i)%q
        end if
      end do
    end if
    
    ! METAR family
    if (associated(y%metar) .and. associated(iv%metar) .and. iv%info(metar)%nlocal > 0) then
      do i = 1, min(size(y%metar), iv%info(metar)%nlocal)
        if (iv%metar(i)%u%qc == 0 .and. abs(iv%metar(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%metar(i)%u
        end if
        if (iv%metar(i)%v%qc == 0 .and. abs(iv%metar(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%metar(i)%v
        end if
        if (iv%metar(i)%t%qc == 0 .and. abs(iv%metar(i)%t%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%metar(i)%t
        end if
        if (iv%metar(i)%p%qc == 0 .and. abs(iv%metar(i)%p%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%metar(i)%p
        end if
        if (iv%metar(i)%q%qc == 0 .and. abs(iv%metar(i)%q%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%metar(i)%q
        end if
      end do
    end if
    
    ! SHIPS family
    if (associated(y%ships) .and. associated(iv%ships) .and. iv%info(ships)%nlocal > 0) then
      do i = 1, min(size(y%ships), iv%info(ships)%nlocal)
        if (iv%ships(i)%u%qc == 0 .and. abs(iv%ships(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ships(i)%u
        end if
        if (iv%ships(i)%v%qc == 0 .and. abs(iv%ships(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ships(i)%v
        end if
        if (iv%ships(i)%t%qc == 0 .and. abs(iv%ships(i)%t%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ships(i)%t
        end if
        if (iv%ships(i)%p%qc == 0 .and. abs(iv%ships(i)%p%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ships(i)%p
        end if
        if (iv%ships(i)%q%qc == 0 .and. abs(iv%ships(i)%q%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ships(i)%q
        end if
      end do
    end if
    
    ! SOUND family - if available and allocated
    ! Note: sound observations have multiple levels per station
    if (associated(y%sound) .and. associated(iv%sound) .and. iv%info(sound)%nlocal > 0) then
      do i = 1, min(size(y%sound), iv%info(sound)%nlocal)
        ! Loop through levels for this station
        do k = 1, iv%info(sound)%levels(i)
          ! Extract U component if QC is good and innovation is non-zero
          if (iv%sound(i)%u(k)%qc == 0 .and. abs(iv%sound(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%sound(i)%u(k)
          end if
          
          ! Extract V component if QC is good and innovation is non-zero
          if (iv%sound(i)%v(k)%qc == 0 .and. abs(iv%sound(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%sound(i)%v(k)
          end if
          
          ! Extract T component if QC is good and innovation is non-zero
          if (iv%sound(i)%t(k)%qc == 0 .and. abs(iv%sound(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%sound(i)%t(k)
          end if
          
          ! Extract Q component if QC is good and innovation is non-zero
          if (iv%sound(i)%q(k)%qc == 0 .and. abs(iv%sound(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%sound(i)%q(k)
          end if
        end do
      end do
    end if
    
    ! BUOY family
    if (associated(y%buoy) .and. associated(iv%buoy) .and. iv%info(buoy)%nlocal > 0) then
      do i = 1, min(size(y%buoy), iv%info(buoy)%nlocal)
        if (iv%buoy(i)%u%qc == 0 .and. abs(iv%buoy(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%buoy(i)%u
        end if
        if (iv%buoy(i)%v%qc == 0 .and. abs(iv%buoy(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%buoy(i)%v
        end if
        if (iv%buoy(i)%t%qc == 0 .and. abs(iv%buoy(i)%t%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%buoy(i)%t
        end if
        if (iv%buoy(i)%p%qc == 0 .and. abs(iv%buoy(i)%p%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%buoy(i)%p
        end if
        if (iv%buoy(i)%q%qc == 0 .and. abs(iv%buoy(i)%q%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%buoy(i)%q
        end if
      end do
    end if
    
    ! AIREP family - multi-level
    if (associated(y%airep) .and. associated(iv%airep) .and. iv%info(airep)%nlocal > 0) then
      do i = 1, min(size(y%airep), iv%info(airep)%nlocal)
        do k = 1, iv%info(airep)%levels(i)
          if (iv%airep(i)%u(k)%qc == 0 .and. abs(iv%airep(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%airep(i)%u(k)
          end if
          if (iv%airep(i)%v(k)%qc == 0 .and. abs(iv%airep(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%airep(i)%v(k)
          end if
          if (iv%airep(i)%t(k)%qc == 0 .and. abs(iv%airep(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%airep(i)%t(k)
          end if
        end do
      end do
    end if
    
    ! PILOT family - multi-level
    if (associated(y%pilot) .and. associated(iv%pilot) .and. iv%info(pilot)%nlocal > 0) then
      do i = 1, min(size(y%pilot), iv%info(pilot)%nlocal)
        do k = 1, iv%info(pilot)%levels(i)
          if (iv%pilot(i)%u(k)%qc == 0 .and. abs(iv%pilot(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%pilot(i)%u(k)
          end if
          if (iv%pilot(i)%v(k)%qc == 0 .and. abs(iv%pilot(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%pilot(i)%v(k)
          end if
        end do
      end do
    end if
    
    ! SONDE_SFC family
    if (associated(y%sonde_sfc) .and. associated(iv%sonde_sfc) .and. iv%info(sonde_sfc)%nlocal > 0) then
      do i = 1, min(size(y%sonde_sfc), iv%info(sonde_sfc)%nlocal)
        if (iv%sonde_sfc(i)%u%qc == 0 .and. abs(iv%sonde_sfc(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%sonde_sfc(i)%u
        end if
        if (iv%sonde_sfc(i)%v%qc == 0 .and. abs(iv%sonde_sfc(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%sonde_sfc(i)%v
        end if
        if (iv%sonde_sfc(i)%t%qc == 0 .and. abs(iv%sonde_sfc(i)%t%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%sonde_sfc(i)%t
        end if
        if (iv%sonde_sfc(i)%p%qc == 0 .and. abs(iv%sonde_sfc(i)%p%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%sonde_sfc(i)%p
        end if
        if (iv%sonde_sfc(i)%q%qc == 0 .and. abs(iv%sonde_sfc(i)%q%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%sonde_sfc(i)%q
        end if
      end do
    end if
    
    ! GEOAMV family - multi-level
    if (associated(y%geoamv) .and. associated(iv%geoamv) .and. iv%info(geoamv)%nlocal > 0) then
      do i = 1, min(size(y%geoamv), iv%info(geoamv)%nlocal)
        do k = 1, iv%info(geoamv)%levels(i)
          if (iv%geoamv(i)%u(k)%qc == 0 .and. abs(iv%geoamv(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%geoamv(i)%u(k)
          end if
          if (iv%geoamv(i)%v(k)%qc == 0 .and. abs(iv%geoamv(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%geoamv(i)%v(k)
          end if
        end do
      end do
    end if
    
    ! POLARAMV family - multi-level
    if (associated(y%polaramv) .and. associated(iv%polaramv) .and. iv%info(polaramv)%nlocal > 0) then
      do i = 1, min(size(y%polaramv), iv%info(polaramv)%nlocal)
        do k = 1, iv%info(polaramv)%levels(i)
          if (iv%polaramv(i)%u(k)%qc == 0 .and. abs(iv%polaramv(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%polaramv(i)%u(k)
          end if
          if (iv%polaramv(i)%v(k)%qc == 0 .and. abs(iv%polaramv(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%polaramv(i)%v(k)
          end if
        end do
      end do
    end if
    
    ! GPSPW family
    if (associated(y%gpspw) .and. associated(iv%gpspw) .and. iv%info(gpspw)%nlocal > 0) then
      do i = 1, min(size(y%gpspw), iv%info(gpspw)%nlocal)
        if (iv%gpspw(i)%tpw%qc == 0 .and. abs(iv%gpspw(i)%tpw%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%gpspw(i)%tpw
        end if
      end do
    end if
    
    ! GPSREF family - multi-level
    if (associated(y%gpsref) .and. associated(iv%gpsref) .and. iv%info(gpsref)%nlocal > 0) then
      do i = 1, min(size(y%gpsref), iv%info(gpsref)%nlocal)
        do k = 1, iv%info(gpsref)%levels(i)
          if (iv%gpsref(i)%ref(k)%qc == 0 .and. abs(iv%gpsref(i)%ref(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%gpsref(i)%ref(k)
          end if
          if (iv%gpsref(i)%p(k)%qc == 0 .and. abs(iv%gpsref(i)%p(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%gpsref(i)%p(k)
          end if
          if (iv%gpsref(i)%t(k)%qc == 0 .and. abs(iv%gpsref(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%gpsref(i)%t(k)
          end if
          if (iv%gpsref(i)%q(k)%qc == 0 .and. abs(iv%gpsref(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%gpsref(i)%q(k)
          end if
        end do
      end do
    end if
    
    ! GPSEPH family - multi-level
    if (associated(y%gpseph) .and. associated(iv%gpseph) .and. iv%info(gpseph)%nlocal > 0) then
      do i = 1, min(size(y%gpseph), iv%info(gpseph)%nlocal)
        do k = 1, iv%info(gpseph)%levels(i)
          if (iv%gpseph(i)%eph(k)%qc == 0 .and. abs(iv%gpseph(i)%eph(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%gpseph(i)%eph(k)
          end if
        end do
      end do
    end if
    
    ! QSCAT family
    if (associated(y%qscat) .and. associated(iv%qscat) .and. iv%info(qscat)%nlocal > 0) then
      do i = 1, min(size(y%qscat), iv%info(qscat)%nlocal)
        if (iv%qscat(i)%u%qc == 0 .and. abs(iv%qscat(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%qscat(i)%u
        end if
        if (iv%qscat(i)%v%qc == 0 .and. abs(iv%qscat(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%qscat(i)%v
        end if
      end do
    end if
    
    ! PROFILER family - multi-level
    if (associated(y%profiler) .and. associated(iv%profiler) .and. iv%info(profiler)%nlocal > 0) then
      do i = 1, min(size(y%profiler), iv%info(profiler)%nlocal)
        do k = 1, iv%info(profiler)%levels(i)
          if (iv%profiler(i)%u(k)%qc == 0 .and. abs(iv%profiler(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%profiler(i)%u(k)
          end if
          if (iv%profiler(i)%v(k)%qc == 0 .and. abs(iv%profiler(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%profiler(i)%v(k)
          end if
        end do
      end do
    end if
    
    ! SSMI_RV family
    if (associated(y%ssmi_rv) .and. associated(iv%ssmi_rv) .and. iv%info(ssmi_rv)%nlocal > 0) then
      do i = 1, min(size(y%ssmi_rv), iv%info(ssmi_rv)%nlocal)
        if (iv%ssmi_rv(i)%Speed%qc == 0 .and. abs(iv%ssmi_rv(i)%Speed%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ssmi_rv(i)%Speed
        end if
        if (iv%ssmi_rv(i)%tpw%qc == 0 .and. abs(iv%ssmi_rv(i)%tpw%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ssmi_rv(i)%tpw
        end if
      end do
    end if
    
    ! SSMI_TB family
    if (associated(y%ssmi_tb) .and. associated(iv%ssmi_tb) .and. iv%info(ssmi_tb)%nlocal > 0) then
      do i = 1, min(size(y%ssmi_tb), iv%info(ssmi_tb)%nlocal)
        if (iv%ssmi_tb(i)%tb19h%qc == 0 .and. abs(iv%ssmi_tb(i)%tb19h%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ssmi_tb(i)%tb19h
        end if
        if (iv%ssmi_tb(i)%tb19v%qc == 0 .and. abs(iv%ssmi_tb(i)%tb19v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ssmi_tb(i)%tb19v
        end if
        if (iv%ssmi_tb(i)%tb22v%qc == 0 .and. abs(iv%ssmi_tb(i)%tb22v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ssmi_tb(i)%tb22v
        end if
        if (iv%ssmi_tb(i)%tb37h%qc == 0 .and. abs(iv%ssmi_tb(i)%tb37h%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ssmi_tb(i)%tb37h
        end if
        if (iv%ssmi_tb(i)%tb37v%qc == 0 .and. abs(iv%ssmi_tb(i)%tb37v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ssmi_tb(i)%tb37v
        end if
        if (iv%ssmi_tb(i)%tb85h%qc == 0 .and. abs(iv%ssmi_tb(i)%tb85h%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ssmi_tb(i)%tb85h
        end if
        if (iv%ssmi_tb(i)%tb85v%qc == 0 .and. abs(iv%ssmi_tb(i)%tb85v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%ssmi_tb(i)%tb85v
        end if
      end do
    end if
    
    ! SSMT1 family - multi-level
    if (associated(y%ssmt1) .and. associated(iv%ssmt1) .and. iv%info(ssmt1)%nlocal > 0) then
      do i = 1, min(size(y%ssmt1), iv%info(ssmt1)%nlocal)
        do k = 1, iv%info(ssmt1)%levels(i)
          if (iv%ssmt1(i)%t(k)%qc == 0 .and. abs(iv%ssmt1(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%ssmt1(i)%t(k)
          end if
        end do
      end do
    end if
    
    ! SSMT2 family - multi-level
    if (associated(y%ssmt2) .and. associated(iv%ssmt2) .and. iv%info(ssmt2)%nlocal > 0) then
      do i = 1, min(size(y%ssmt2), iv%info(ssmt2)%nlocal)
        do k = 1, iv%info(ssmt2)%levels(i)
          if (iv%ssmt2(i)%rh(k)%qc == 0 .and. abs(iv%ssmt2(i)%rh(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%ssmt2(i)%rh(k)
          end if
        end do
      end do
    end if
    
    ! SATEM family - multi-level
    if (associated(y%satem) .and. associated(iv%satem) .and. iv%info(satem)%nlocal > 0) then
      do i = 1, min(size(y%satem), iv%info(satem)%nlocal)
        do k = 1, iv%info(satem)%levels(i)
          if (iv%satem(i)%thickness(k)%qc == 0 .and. abs(iv%satem(i)%thickness(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%satem(i)%thickness(k)
          end if
        end do
      end do
    end if
    
    ! PSEUDO family
    if (associated(y%pseudo) .and. associated(iv%pseudo) .and. iv%info(pseudo)%nlocal > 0) then
      do i = 1, min(size(y%pseudo), iv%info(pseudo)%nlocal)
        if (iv%pseudo(i)%u%qc == 0 .and. abs(iv%pseudo(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%pseudo(i)%u
        end if
        if (iv%pseudo(i)%v%qc == 0 .and. abs(iv%pseudo(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%pseudo(i)%v
        end if
        if (iv%pseudo(i)%t%qc == 0 .and. abs(iv%pseudo(i)%t%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%pseudo(i)%t
        end if
        if (iv%pseudo(i)%p%qc == 0 .and. abs(iv%pseudo(i)%p%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%pseudo(i)%p
        end if
        if (iv%pseudo(i)%q%qc == 0 .and. abs(iv%pseudo(i)%q%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%pseudo(i)%q
        end if
      end do
    end if
    
    ! BOGUS family - multi-level
    if (associated(y%bogus) .and. associated(iv%bogus) .and. iv%info(bogus)%nlocal > 0) then
      do i = 1, min(size(y%bogus), iv%info(bogus)%nlocal)
        do k = 1, iv%info(bogus)%levels(i)
          if (iv%bogus(i)%u(k)%qc == 0 .and. abs(iv%bogus(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%bogus(i)%u(k)
          end if
          if (iv%bogus(i)%v(k)%qc == 0 .and. abs(iv%bogus(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%bogus(i)%v(k)
          end if
          if (iv%bogus(i)%t(k)%qc == 0 .and. abs(iv%bogus(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%bogus(i)%t(k)
          end if
          if (iv%bogus(i)%q(k)%qc == 0 .and. abs(iv%bogus(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%bogus(i)%q(k)
          end if
        end do
        if (iv%bogus(i)%slp%qc == 0 .and. abs(iv%bogus(i)%slp%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%bogus(i)%slp
        end if
      end do
    end if
    
    ! AIRSR family - multi-level
    if (associated(y%airsr) .and. associated(iv%airsr) .and. iv%info(airsr)%nlocal > 0) then
      do i = 1, min(size(y%airsr), iv%info(airsr)%nlocal)
        do k = 1, iv%info(airsr)%levels(i)
          if (iv%airsr(i)%t(k)%qc == 0 .and. abs(iv%airsr(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%airsr(i)%t(k)
          end if
          if (iv%airsr(i)%q(k)%qc == 0 .and. abs(iv%airsr(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%airsr(i)%q(k)
          end if
        end do
      end do
    end if
    
    ! MTGIRS family - multi-level
    if (associated(y%mtgirs) .and. associated(iv%mtgirs) .and. iv%info(mtgirs)%nlocal > 0) then
      do i = 1, min(size(y%mtgirs), iv%info(mtgirs)%nlocal)
        do k = 1, iv%info(mtgirs)%levels(i)
          if (iv%mtgirs(i)%u(k)%qc == 0 .and. abs(iv%mtgirs(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%mtgirs(i)%u(k)
          end if
          if (iv%mtgirs(i)%v(k)%qc == 0 .and. abs(iv%mtgirs(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%mtgirs(i)%v(k)
          end if
          if (iv%mtgirs(i)%t(k)%qc == 0 .and. abs(iv%mtgirs(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%mtgirs(i)%t(k)
          end if
          if (iv%mtgirs(i)%q(k)%qc == 0 .and. abs(iv%mtgirs(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%mtgirs(i)%q(k)
          end if
        end do
      end do
    end if
    
    ! TAMDAR family - multi-level
    if (associated(y%tamdar) .and. associated(iv%tamdar) .and. iv%info(tamdar)%nlocal > 0) then
      do i = 1, min(size(y%tamdar), iv%info(tamdar)%nlocal)
        do k = 1, iv%info(tamdar)%levels(i)
          if (iv%tamdar(i)%u(k)%qc == 0 .and. abs(iv%tamdar(i)%u(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%tamdar(i)%u(k)
          end if
          if (iv%tamdar(i)%v(k)%qc == 0 .and. abs(iv%tamdar(i)%v(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%tamdar(i)%v(k)
          end if
          if (iv%tamdar(i)%t(k)%qc == 0 .and. abs(iv%tamdar(i)%t(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%tamdar(i)%t(k)
          end if
          if (iv%tamdar(i)%q(k)%qc == 0 .and. abs(iv%tamdar(i)%q(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%tamdar(i)%q(k)
          end if
        end do
      end do
    end if
    
    ! TAMDAR_SFC family
    if (associated(y%tamdar_sfc) .and. associated(iv%tamdar_sfc) .and. iv%info(tamdar_sfc)%nlocal > 0) then
      do i = 1, min(size(y%tamdar_sfc), iv%info(tamdar_sfc)%nlocal)
        if (iv%tamdar_sfc(i)%u%qc == 0 .and. abs(iv%tamdar_sfc(i)%u%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%tamdar_sfc(i)%u
        end if
        if (iv%tamdar_sfc(i)%v%qc == 0 .and. abs(iv%tamdar_sfc(i)%v%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%tamdar_sfc(i)%v
        end if
        if (iv%tamdar_sfc(i)%t%qc == 0 .and. abs(iv%tamdar_sfc(i)%t%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%tamdar_sfc(i)%t
        end if
        if (iv%tamdar_sfc(i)%p%qc == 0 .and. abs(iv%tamdar_sfc(i)%p%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%tamdar_sfc(i)%p
        end if
        if (iv%tamdar_sfc(i)%q%qc == 0 .and. abs(iv%tamdar_sfc(i)%q%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%tamdar_sfc(i)%q
        end if
      end do
    end if
    
    ! RAIN family
    if (associated(y%rain) .and. associated(iv%rain) .and. iv%info(rain)%nlocal > 0) then
      do i = 1, min(size(y%rain), iv%info(rain)%nlocal)
        if (iv%rain(i)%rain%qc == 0 .and. abs(iv%rain(i)%rain%inv) > 1.0e-10) then
          count = count + 1
          observations(count) = y%rain(i)%rain
        end if
      end do
    end if
    
    ! RADAR family - multi-level
    if (associated(y%radar) .and. associated(iv%radar) .and. iv%info(radar)%nlocal > 0) then
      do i = 1, min(size(y%radar), iv%info(radar)%nlocal)
        do k = 1, iv%info(radar)%levels(i)
          if (iv%radar(i)%rv(k)%qc == 0 .and. abs(iv%radar(i)%rv(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%radar(i)%rv(k)
          end if
          if (iv%radar(i)%rf(k)%qc == 0 .and. abs(iv%radar(i)%rf(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%radar(i)%rf(k)
          end if
        end do
      end do
    end if
    
    ! LIGHTNING family - multi-level
    if (associated(y%lightning) .and. associated(iv%lightning) .and. iv%info(lightning)%nlocal > 0) then
      do i = 1, min(size(y%lightning), iv%info(lightning)%nlocal)
        do k = 1, iv%info(lightning)%levels(i)
          if (iv%lightning(i)%w(k)%qc == 0 .and. abs(iv%lightning(i)%w(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%lightning(i)%w(k)
          end if
          if (iv%lightning(i)%div(k)%qc == 0 .and. abs(iv%lightning(i)%div(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%lightning(i)%div(k)
          end if
          if (iv%lightning(i)%qv(k)%qc == 0 .and. abs(iv%lightning(i)%qv(k)%inv) > 1.0e-10) then
            count = count + 1
            observations(count) = y%lightning(i)%qv(k)
          end if
        end do
      end do
    end if
    
    num_observations = count
    wrfda_extract_observations = 0
    
  end function wrfda_extract_observations


  ! Initialize WRFDA variables for 3D-Var analysis
  !> @brief Copy namelist configuration from model_config_rec to da_control module
  !> 
  !> @details This function extracts and reuses the config copy logic from WRFDA's
  !>          da_wrfvar_init1 subroutine (lines 87-92 of da_wrfvar_init1.inc).
  !>          
  !>          We cannot call da_wrfvar_init1 directly because it also:
  !>          - Calls init_modules (already done in WRFConfigManager)
  !>          - Calls wrfu_initialize (already done)
  !>          - Calls initial_config to read namelist (already done)
  !>          - Sets up MPI variables (already done)
  !>          
  !>          This function extracts ONLY the config copy step using the EXACT
  !>          same mechanism as da_wrfvar_init1: including config_assigns.inc
  !>          with SOURCE_RECORD defined as model_config_rec%.
  !>          
  !> @note This must be called after initial_config() to populate da_control variables
  !> @see da_wrfvar_init1.inc lines 87-92 in WRFDA source
  subroutine copy_config_to_da_control() bind(C, name="copy_config_to_da_control")
    use module_configure, only: model_config_rec
    use da_control
    implicit none
    
    ! Copy namelist variables from model_config_rec to da_control module
    ! This is the EXACT code from da_wrfvar_init1.inc lines 87-92:
    !   #define SOURCE_RECORD model_config_rec%
    !   #define DEST_RECORD
    !   #include "config_assigns.inc"
    ! 
    ! config_assigns.inc is generated by WRF's Registry and contains assignments like:
    !   use_synopobs = model_config_rec%use_synopobs
    !   use_metarobs = model_config_rec%use_metarobs
    !   ... (2000+ lines of similar assignments)
    
#define SOURCE_RECORD model_config_rec%
#define DEST_RECORD
#include "config_assigns.inc"
#undef SOURCE_RECORD
#undef DEST_RECORD
    
  end subroutine copy_config_to_da_control

  !> @brief Validate WRFDA configuration for common conflicts and errors
  !> 
  !> @details This function performs the same sanity checks as WRFDA's da_solve.inc
  !>          to catch configuration errors early. These checks replicate the logic
  !>          from da_solve.inc lines 168-252 (Initial checks section).
  !>          
  !> @note This must be called after copy_config_to_da_control()
  !> @see da_solve.inc lines 168-252 in WRFDA source
  !> @return Integer error code (0 = success, non-zero = validation error)
  subroutine validate_wrfda_config(error_code) bind(C, name="validate_wrfda_config")
    use da_control, only: use_gpsrefobs, use_gpsephobs, use_radar_rf, use_radar_rhv, &
                          use_radarobs, cv_options_hum, cv_options_hum_specific_humidity, &
                          cv_options_hum_relative_humidity, vert_corr, vert_corr_2, &
                          vertical_ip, vertical_ip_0, vertical_ip_delta_p, &
                          cv_options, cloud_cv_options, ensdim_alpha, alphacv_method, &
                          alphacv_method_xa, anal_type_hybrid_dual_res, radar_rf_opt, &
                          alpha_hydrometeors
    use da_wrf_interfaces, only: wrf_message, wrf_error_fatal
    implicit none
    
    integer(c_int), intent(out) :: error_code
    character(len=512) :: message
    
    error_code = 0_c_int
    
    ! ========================================================================
    ! GPS Observation Conflicts (da_solve.inc lines 248-252)
    ! ========================================================================
    if (use_gpsrefobs .and. use_gpsephobs) then
      call wrf_message("ERROR: Configuration Validation Failed")
      call wrf_message("Both use_gpsrefobs and use_gpsephobs are set to true")
      call wrf_message("You must choose EITHER use_gpsrefobs OR use_gpsephobs, not both")
      call wrf_message("Please update your namelist.input &wrfvar4 section and rerun")
      error_code = 1_c_int
      return
    end if
    
    ! ========================================================================
    ! Radar Observation Conflicts (da_solve.inc lines 225-238)
    ! ========================================================================
    if (use_radarobs) then
      if (use_radar_rf .and. use_radar_rhv) then
        call wrf_message("ERROR: Configuration Validation Failed")
        call wrf_message("Both use_radar_rf and use_radar_rhv are set to true")
        call wrf_message("You must choose EITHER use_radar_rf OR use_radar_rhv, not both")
        call wrf_message("Please update your namelist.input &wrfvar4 section and rerun")
        error_code = 2_c_int
        return
      end if
      
      if (use_radar_rf .and. radar_rf_opt == 1 .and. cloud_cv_options /= 1) then
        call wrf_message("ERROR: Configuration Validation Failed")
        call wrf_message("For use_radar_rf with radar_rf_opt=1, you must set cloud_cv_options=1")
        call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
        error_code = 3_c_int
        return
      end if
      
      if (use_radar_rhv .and. cloud_cv_options == 1) then
        call wrf_message("ERROR: Configuration Validation Failed")
        call wrf_message("For use_radar_rhv, you must set cloud_cv_options=2 or 3")
        call wrf_message("(cloud_cv_options=2 requires cloudy be.dat)")
        call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
        error_code = 4_c_int
        return
      end if
    end if
    
    ! ========================================================================
    ! Control Variable Options Validation (da_solve.inc lines 168-173)
    ! ========================================================================
    if (cv_options_hum /= cv_options_hum_specific_humidity .and. &
        cv_options_hum /= cv_options_hum_relative_humidity) then
      write(message, '(A,I3)') &
        "ERROR: Invalid cv_options_hum = ", cv_options_hum
      call wrf_message(trim(message))
      call wrf_message("cv_options_hum must be 1 (specific humidity) or 2 (relative humidity)")
      call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
      error_code = 5_c_int
      return
    end if
    
    ! ========================================================================
    ! Vertical Interpolation Validation (da_solve.inc lines 175-181)
    ! ========================================================================
    if (vert_corr == vert_corr_2) then
      if (vertical_ip < vertical_ip_0 .or. vertical_ip > vertical_ip_delta_p) then
        write(message, '(A,I3)') &
          "ERROR: Invalid vertical_ip = ", vertical_ip
        call wrf_message(trim(message))
        write(message, '(A,I1,A,I1)') &
          "vertical_ip must be between ", vertical_ip_0, " and ", vertical_ip_delta_p
        call wrf_message(trim(message))
        call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
        error_code = 6_c_int
        return
      end if
    end if
    
    ! ========================================================================
    ! CV Options = 3 Constraints (da_solve.inc lines 206-221)
    ! ========================================================================
    if (cv_options == 3) then
      if (ensdim_alpha > 0) then
        call wrf_message("ERROR: Configuration Validation Failed")
        call wrf_message("Alpha control variables are not implemented for cv_options=3")
        call wrf_message("Please set cv_options=5 or 7, OR set ensdim_alpha=0")
        call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
        error_code = 7_c_int
        return
      end if
    end if
    
    ! ========================================================================
    ! Dual-Resolution Hybrid Constraints (da_solve.inc lines 194-198)
    ! ========================================================================
    if (anal_type_hybrid_dual_res .and. alphacv_method /= alphacv_method_xa) then
      call wrf_message("ERROR: Configuration Validation Failed")
      call wrf_message("Dual-resolution hybrid requires alphacv_method=2")
      write(message, '(A,I1)') "Current alphacv_method = ", alphacv_method
      call wrf_message(trim(message))
      call wrf_message("Please update your namelist.input &wrfvar17 section and rerun")
      error_code = 8_c_int
      return
    end if
    
    ! ========================================================================
    ! Alpha Hydrometeors Constraints (da_solve.inc lines 240-245)
    ! ========================================================================
    if (ensdim_alpha > 0 .and. alpha_hydrometeors) then
      if (cloud_cv_options == 1) then
        call wrf_message("ERROR: Configuration Validation Failed")
        call wrf_message("alpha_hydrometeors is not implemented for cloud_cv_options=1")
        call wrf_message("Please set cloud_cv_options=3")
        call wrf_message("Please update your namelist.input &wrfvar7 section and rerun")
        error_code = 9_c_int
        return
      end if
    end if
    
    ! All checks passed
    call wrf_message("Configuration validation passed - all settings are consistent")
    
  end subroutine validate_wrfda_config

  subroutine initialize_wrfda_3dvar() bind(C, name="initialize_wrfda_3dvar")
    implicit none
    
    ! Set variables for 3D-Var analysis
    ! var4d_run = .false. indicates 3D-Var (not 4D-Var)
    ! num_fgat_time = 1 indicates single time slot for 3D-Var
    ! sfc_assi_options = sfc_assi_options_1 sets surface assimilation options
    var4d_run = .false.
    num_fgat_time = 1
    sfc_assi_options = sfc_assi_options_1
    
  end subroutine initialize_wrfda_3dvar
  

  ! Initialize domain clock for WRFDA time management
  subroutine initialize_domain_clock(grid)
    implicit none
    type(domain), intent(inout) :: grid
    character(len=256) :: timestr
    logical, save :: wrfu_initialized = .false.
    type(WRFU_Time) :: start_time, stop_time
    type(WRFU_TimeInterval) :: time_step
    integer :: rc
    
    ! Initialize WRFU (required for domain_clock_set calls)
    if (.not. wrfu_initialized) then
      call wrfu_initialize(defaultCalKind=wrfu_cal_gregorian)
      wrfu_initialized = .true.
    end if
    
    ! Set a default analysis time for the domain clock
    timestr = '2000-01-24_12:00:00'  ! Default analysis time
    
    ! Convert time string to WRFU_Time using WRF's wrf_atotime function
    call wrf_atotime(timestr, start_time)
    call wrf_atotime(timestr, stop_time)
    
    ! Create time step interval (0 seconds for analysis)
    call WRFU_TimeIntervalSet(time_step, S=0, rc=rc)
    if (rc /= WRFU_SUCCESS) then
      print *, "WRFDA ERROR: Failed to create time step interval"
      return
    end if
    
    ! Create the domain clock using WRFU_ClockCreate
    ! Note: domain_clock is a TYPE(WRFU_Clock), not a pointer, so no allocation needed
    grid%domain_clock = WRFU_ClockCreate(TimeStep=time_step, StartTime=start_time, StopTime=stop_time, rc=rc)
    if (rc /= WRFU_SUCCESS) then
      print *, "WRFDA ERROR: Failed to create domain clock"
      return
    end if
    
    ! Set the domain_clock_created flag to true
    grid%domain_clock_created = .true.
    
  end subroutine initialize_domain_clock

  ! =============================================================================
  ! NEW: WRFDA First Guess Initialization Functions
  ! These functions implement the refactored WRFState initialization strategy
  ! that leverages WRFDA's proven initialization pipeline
  ! =============================================================================

  !> @brief Load WRF first guess using WRFDA's da_med_initialdata_input and da_setup_firstguess
  !> @details This function implements the complete WRFDA initialization pipeline:
  !>          1. Calls da_med_initialdata_input to read NetCDF file into grid structure
  !>          2. Calls da_setup_firstguess_wrf to setup grid and call da_transfer_wrftoxb
  !>          3. All WRFDA diagnostics and derived fields are computed properly
  !> @param[in] filename C string containing path to WRF NetCDF file
  !> @param[in] filename_len Length of filename string
  !> @return Integer status code (0 = success, non-zero = error)
!> @brief Update da_control module variables from grid structure
!> @details Sets all grid-dependent da_control module variables after da_med_initialdata_input.
!>          This includes:
!>            (1) Vertical coordinate coefficients (c1f, c2f, c3f, c4f, c1h, c2h, c3h, c4h)
!>            (2) Base state parameters (base_pres, base_temp, base_lapse, iso_temp, etc.)
!>          Map projection variables are set separately by da_setup_firstguess_wrf.
!>          This function must be called AFTER da_med_initialdata_input and BEFORE
!>          da_setup_firstguess_wrf/da_transfer_wrftoxb.
!> @param[in] grid Pointer to WRFDA grid structure with populated state from NetCDF
subroutine update_da_control_from_grid(grid)
  use module_domain, only: domain
  use da_control, only: c1f, c2f, c3f, c4f, c1h, c2h, c3h, c4h, &
                        base_pres, base_temp, base_lapse, iso_temp, &
                        base_pres_strat, base_lapse_strat
  use da_wrf_interfaces, only: wrf_message
  implicit none
  
  type(domain), pointer, intent(in) :: grid
  integer :: kms, kme
  character(len=256) :: msg
  
  call wrf_message("======== Updating da_control module variables from grid ========")
  
  ! -------------------------------------------------------------------------
  ! [1] Update vertical coordinate coefficients
  ! -------------------------------------------------------------------------
  kms = grid%sm33
  kme = grid%em33
  
  ! Allocate da_control module arrays if not already allocated
  if (.not. allocated(c1f)) allocate(c1f(kms:kme))
  if (.not. allocated(c2f)) allocate(c2f(kms:kme))
  if (.not. allocated(c3f)) allocate(c3f(kms:kme))
  if (.not. allocated(c4f)) allocate(c4f(kms:kme))
  if (.not. allocated(c1h)) allocate(c1h(kms:kme))
  if (.not. allocated(c2h)) allocate(c2h(kms:kme))
  if (.not. allocated(c3h)) allocate(c3h(kms:kme))
  if (.not. allocated(c4h)) allocate(c4h(kms:kme))
  
  write(msg, '(A,I3,A,I3,A)') "[1] Vertical coordinates (kms:kme = ", kms, ":", kme, ")"
  call wrf_message(msg)
  
  ! Copy from grid structure (populated by da_med_initialdata_input)
  c1f = grid%c1f
  c2f = grid%c2f
  c3f = grid%c3f
  c4f = grid%c4f
  c1h = grid%c1h
  c2h = grid%c2h
  c3h = grid%c3h
  c4h = grid%c4h
  
  ! Handle non-hybrid coordinates (backward compatibility)
  ! For input files prior to V3.9, grid%hybrid_opt is set to 0 by da_med_initialdata_input
  if (grid%hybrid_opt <= 0) then
    write(msg, '(A,I2,A)') "    Hybrid_opt = ", grid%hybrid_opt, " -> using pure eta coordinates"
    call wrf_message(msg)
    
    ! Fall back to pure eta coordinates
    c3f = grid%znw  ! Eta levels on full (w) layers
    c3h = grid%znu  ! Eta levels on half (mass) layers
    c4f = 0.0
    c4h = 0.0
    c1f = 1.0
    c1h = 1.0
    c2f = 0.0
    c2h = 0.0
  else
    write(msg, '(A,I2,A)') "    Hybrid_opt = ", grid%hybrid_opt, " -> using hybrid coordinates"
    call wrf_message(msg)
  end if
  
  ! -------------------------------------------------------------------------
  ! [2] Update base state parameters
  ! -------------------------------------------------------------------------
  write(msg, '(A)') "[2] Base state parameters"
  call wrf_message(msg)
  
  base_pres  = grid%p00
  base_temp  = grid%t00
  base_lapse = grid%tlp
  iso_temp   = grid%tiso
  base_pres_strat  = grid%p_strat
  base_lapse_strat = grid%tlp_strat
  
  write(msg, '(A,F10.1,A)') "    base_pres        = ", base_pres, " Pa"
  call wrf_message(msg)
  write(msg, '(A,F10.2,A)') "    base_temp        = ", base_temp, " K"
  call wrf_message(msg)
  write(msg, '(A,F10.6)') "    base_lapse       = ", base_lapse
  call wrf_message(msg)
  write(msg, '(A,F10.2,A)') "    iso_temp         = ", iso_temp, " K"
  call wrf_message(msg)
  write(msg, '(A,F10.1,A)') "    base_pres_strat  = ", base_pres_strat, " Pa"
  call wrf_message(msg)
  write(msg, '(A,F10.6)') "    base_lapse_strat = ", base_lapse_strat
  call wrf_message(msg)
  
  ! Validate - base_temp should be > 100 K and base_pres > 10000 Pa
  if ( base_temp < 100.0 .or. base_pres < 10000.0 ) then
    write(msg, '(A)') "ERROR: Base state parameters not found in NetCDF file!"
    call wrf_message(msg)
    write(msg, '(A)') "Add use_baseparam_fr_nml = .true. in namelist.input &dynamics"
    call wrf_message(msg)
    write(msg, '(A,F10.2,A,F10.1)') "Got: base_temp = ", base_temp, " K, base_pres = ", base_pres
    call wrf_message(msg)
    stop "FATAL: Missing or invalid base state parameters"
  end if
  
  call wrf_message("======== da_control module variables updated successfully ========")
  
end subroutine update_da_control_from_grid

!> @brief Cleanup da_control module vertical coordinate variables
!> @details Deallocates vertical coordinate arrays in da_control module.
!>          WARNING: These are module-level variables shared across ALL grid instances!
!>          Should ONLY be called during final WRFDA shutdown, NOT in instance destructors.
!>          Calling this prematurely will break other WRFState/WRFGeometry instances.
subroutine cleanup_da_control_vertical_coords()
  use da_control, only: c1f, c2f, c3f, c4f, c1h, c2h, c3h, c4h
  use da_wrf_interfaces, only: wrf_message
  implicit none
  
  call wrf_message("Cleaning up da_control vertical coordinates")
  
  if (allocated(c1f)) deallocate(c1f)
  if (allocated(c2f)) deallocate(c2f)
  if (allocated(c3f)) deallocate(c3f)
  if (allocated(c4f)) deallocate(c4f)
  if (allocated(c1h)) deallocate(c1h)
  if (allocated(c2h)) deallocate(c2h)
  if (allocated(c3h)) deallocate(c3h)
  if (allocated(c4h)) deallocate(c4h)
  
  call wrf_message("da_control vertical coordinates cleaned up")
  
end subroutine cleanup_da_control_vertical_coords

!> @brief C-callable wrapper for cleanup_da_control_vertical_coords
!> @details Allows C++ code to cleanup da_control module vertical coordinates.
!>          WARNING: Only call during final program shutdown, NOT in destructors!
!>          These are shared module-level variables across all instances.
subroutine wrfda_cleanup_vertical_coords() bind(C, name="wrfda_cleanup_vertical_coords")
  implicit none
  
  call cleanup_da_control_vertical_coords()
  
end subroutine wrfda_cleanup_vertical_coords

integer(c_int) function wrfda_load_first_guess(grid_ptr, filename, filename_len) &
    bind(C, name="wrfda_load_first_guess")
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_char, c_int, c_null_char, c_associated, c_double
  use module_configure, only: grid_config_rec_type, model_config_rec, model_to_grid_config_rec
  use module_domain, only: domain, head_grid
  use module_tiles, only: set_tiles
  use da_wrfvar_io, only: da_med_initialdata_input
  use da_transfer_model, only: da_setup_firstguess_wrf
  use da_wrf_interfaces, only: wrf_message
  implicit none
  
  ! Input parameters
  type(c_ptr), value, intent(in) :: grid_ptr
  type(c_ptr), value, intent(in) :: filename
  integer(c_int), value, intent(in) :: filename_len
  
  ! Local variables
  type(domain), pointer :: grid
  character(len=512) :: fg_filename
  type(xbx_type) :: xbx
  type(grid_config_rec_type) :: config_flags
  integer :: i, max_len
  character(len=256) :: msg
  character(kind=c_char), pointer :: fptr(:)
  
  ! Initialize return code
  wrfda_load_first_guess = 0
  
  ! Validate grid pointer (must be pre-allocated via wrfda_alloc_and_init_domain)
  if (.not. c_associated(grid_ptr)) then
    call wrf_message("ERROR: grid_ptr is null. Allocate domain first via wrfda_alloc_and_init_domain")
    wrfda_load_first_guess = -1
    return
  end if
  
  ! Convert C pointer to Fortran pointer (this is head_grid from WRFConfigManager)
  call c_f_pointer(grid_ptr, grid)
  
  ! Convert C filename string to Fortran string
  if (filename_len < 0 .or. filename_len > 1024) then
    call wrf_message("ERROR: Invalid filename_len")
    wrfda_load_first_guess = -1
    return
  end if
  
  fg_filename = ""
  call c_f_pointer(filename, fptr, [filename_len])
  max_len = min(filename_len, len(fg_filename))
  do i = 1, max_len
    if (fptr(i) == c_null_char) exit
    fg_filename(i:i) = fptr(i)
  end do
  
  write(msg, '(A,A)') "WRFDA: Loading first guess from file: ", trim(fg_filename)
  call wrf_message(msg)
  
  ! Extract domain-specific config in a single canonical call
  ! Note: This was already called in wrfda_alloc_and_init_domain, but we need
  ! a local config_flags for da_med_initialdata_input
  call model_to_grid_config_rec(head_grid%id, model_config_rec, config_flags)
  
  ! Call da_med_initialdata_input to read NetCDF file
  call da_med_initialdata_input(head_grid, config_flags, trim(fg_filename))
  
  ! Update da_control module variables (vertical coords + base state parameters)
  ! This must be done AFTER da_med_initialdata_input (which populates grid structure)
  ! and BEFORE da_setup_firstguess_wrf/da_transfer_wrftoxb (which use these module variables)
  call update_da_control_from_grid(head_grid)
  
  ! Initialize WRFDA module-level variables following standard WRFDA sequence
  ! This matches the exact sequence used in da_solve_init.inc (lines 33-40)
  ! Step 1: De-reference dimension information from grid structure
  !         Sets module variables: ids, ide, jds, jde, ips, ipe, jps, jpe, kms, kme, etc.
  call da_copy_dims(head_grid)
  
  ! Step 2: Compute tile starting/stopping locations
  !         Sets grid%i_start, grid%i_end, grid%j_start, grid%j_end, grid%num_tiles
  call set_tiles(head_grid, ids, ide, jds, jde, ips, ipe, jps, jpe)
  
  ! Step 3: Copy tile dimensions to module-level variables
  !         Sets module variables: its, ite, jts, jte, kts, kte
  call da_copy_tile_dims(head_grid)
  
  ! Set surface assimilation options
  sfc_assi_options = sfc_assi_options_1
  
  ! Call da_setup_firstguess_wrf to setup grid and call da_transfer_wrftoxb
  call da_setup_firstguess_wrf(xbx, head_grid, config_flags, .false.)
  
  write(msg, '(A)') "WRFDA: First guess loaded successfully"
  call wrf_message(msg)
    
  end function wrfda_load_first_guess

  !> @brief Extract background state (xb) data from WRFDA grid to C++ arrays
  !> @details Extracts the fully initialized background state from WRFDA's grid structure
  !>          to flat C++ arrays for use in the WRFState class
  !> @param[out] u U-wind component (3D: nx*ny*nz)
  !> @param[out] v V-wind component (3D: nx*ny*nz)
  !> @param[out] t Temperature perturbation (3D: nx*ny*nz)
  !> @param[out] q Specific humidity (3D: nx*ny*nz)
  !> @param[out] psfc Surface pressure (2D: nx*ny)
  !> @param[out] p Full pressure (3D: nx*ny*nz)
  !> @param[out] ph Geopotential perturbation (3D staggered: nx*ny*(nz+1))
  !> @param[out] phb Base state geopotential (3D staggered: nx*ny*(nz+1))
  !> @param[out] hgt Terrain height (2D: nx*ny)
  !> @param[out] lats 2D latitude array (2D: nx*ny)
  !> @param[out] lons 2D longitude array (2D: nx*ny)
  !> @param[out] nx Grid x-dimension
  !> @param[out] ny Grid y-dimension
  !> @param[out] nz Grid z-dimension
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_extract_background_state( &
      grid_ptr, u, v, t, q, psfc, p, ph, phb, hgt, lats, lons, &
      nx, ny, nz) bind(C, name="wrfda_extract_background_state")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use module_domain, only: domain
    implicit none
    
    ! Input parameters
    type(c_ptr), value, intent(in) :: grid_ptr
    
    ! Output parameters
    real(c_double), intent(out) :: u(*), v(*), t(*), q(*), psfc(*)
    real(c_double), intent(out) :: p(*), ph(*), phb(*), hgt(*)
    real(c_double), intent(out) :: lats(*), lons(*)
    integer(c_int), intent(out) :: nx, ny, nz
    
    ! Local variables
    integer :: i, j, k, idx, idx_3d, staggered_nz
    type(domain), pointer :: grid
    
    ! Initialize return code
    wrfda_extract_background_state = 0
    
    ! Handle grid pointer: use persistent_grid if grid_ptr is null
    if (.not. c_associated(grid_ptr)) then
      if (.not. grid_initialized .or. .not. associated(persistent_grid)) then
        wrfda_extract_background_state = -1
        return
      end if
      grid => persistent_grid
    else
      call c_f_pointer(grid_ptr, grid)
    end if
    
    ! Get grid dimensions
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    staggered_nz = nz + 1
    
    ! Check if xb is allocated
    if (.not. associated(grid%xb%u) .or. .not. associated(grid%xb%v) .or. &
        .not. associated(grid%xb%t) .or. .not. associated(grid%xb%q) .or. &
        .not. associated(grid%xb%psfc)) then
      wrfda_extract_background_state = 2
      return
    end if
    
    ! Extract 3D fields (U, V, T, Q, P)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! C++ row-major indexing: [i][j][k] -> i + j*nx + k*nx*ny
          ! But we're using 1-based Fortran indexing, so: (i-1) + (j-1)*nx + (k-1)*nx*ny
          idx = (i-1) + (j-1)*nx + (k-1)*nx*ny + 1
          u(idx) = real(grid%xb%u(i,j,k), kind=c_double)
          v(idx) = real(grid%xb%v(i,j,k), kind=c_double)
          t(idx) = real(grid%xb%t(i,j,k), kind=c_double)
          q(idx) = real(grid%xb%q(i,j,k), kind=c_double)
          p(idx) = real(grid%xb%p(i,j,k), kind=c_double)
        end do
      end do
    end do
    
    ! Extract 2D surface fields (PSFC, HGT)
    do j = 1, ny
      do i = 1, nx
        idx = (i-1) + (j-1)*nx + 1
        psfc(idx) = real(grid%xb%psfc(i,j), kind=c_double)
        hgt(idx) = real(grid%xb%terr(i,j), kind=c_double)
        lats(idx) = real(grid%xb%lat(i,j), kind=c_double)
        lons(idx) = real(grid%xb%lon(i,j), kind=c_double)
      end do
    end do
    
    ! Extract vertically staggered fields (PH, PHB) if available
    if (associated(grid%ph_2) .and. associated(grid%phb)) then
      do k = 1, staggered_nz
        do j = 1, ny
          do i = 1, nx
            idx_3d = (i-1) + (j-1)*nx + (k-1)*nx*ny + 1
            ph(idx_3d) = real(grid%ph_2(i,j,k), kind=c_double)
            phb(idx_3d) = real(grid%phb(i,j,k), kind=c_double)
          end do
        end do
      end do
    end if
    
  end function wrfda_extract_background_state

  !> @brief Extract additional WRF fields needed for full state representation
  !> @details Extracts additional fields beyond core state variables
  !> @param[out] w Vertical velocity (3D staggered: nx*ny*(nz+1))
  !> @param[out] mu Dry air mass perturbation (2D: nx*ny)
  !> @param[out] mub Base state dry air mass (2D: nx*ny)
  !> @param[out] pb Base state pressure (3D: nx*ny*nz)
  !> @param[out] t_init Initial temperature field (3D: nx*ny*nz)
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_extract_additional_fields( &
      grid_ptr, w, mu, mub, pb, t_init) bind(C, name="wrfda_extract_additional_fields")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use module_domain, only: domain
    implicit none
    
    ! Input parameters
    type(c_ptr), value, intent(in) :: grid_ptr
    
    ! Output parameters
    real(c_double), intent(out) :: w(*), mu(*), mub(*), pb(*), t_init(*)
    
    ! Local variables
    integer :: i, j, k, idx, idx_3d, nx, ny, nz, staggered_nz
    type(domain), pointer :: grid
    
    ! Initialize return code
    wrfda_extract_additional_fields = 0
    
    ! Handle grid pointer: use persistent_grid if grid_ptr is null
    if (.not. c_associated(grid_ptr)) then
      if (.not. grid_initialized .or. .not. associated(persistent_grid)) then
        wrfda_extract_additional_fields = -1
        return
      end if
      grid => persistent_grid
    else
      call c_f_pointer(grid_ptr, grid)
    end if
    
    ! Get grid dimensions
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    staggered_nz = nz + 1
    
    ! Extract W field (vertically staggered)
    if (associated(grid%w_2)) then
      do k = 1, staggered_nz
        do j = 1, ny
          do i = 1, nx
            idx_3d = (i-1) + (j-1)*nx + (k-1)*nx*ny + 1
            w(idx_3d) = real(grid%w_2(i,j,k), kind=c_double)
          end do
        end do
      end do
    end if
    
    ! Extract MU and MUB fields (2D)
    if (associated(grid%mu_2) .and. associated(grid%mub)) then
      do j = 1, ny
        do i = 1, nx
          idx = (i-1) + (j-1)*nx + 1
          mu(idx) = real(grid%mu_2(i,j), kind=c_double)
          mub(idx) = real(grid%mub(i,j), kind=c_double)
        end do
      end do
    end if
    
    ! Extract PB and T_INIT fields (3D)
    if (associated(grid%pb) .and. associated(grid%t_init)) then
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            idx = (i-1) + (j-1)*nx + (k-1)*nx*ny + 1
            pb(idx) = real(grid%pb(i,j,k), kind=c_double)
            t_init(idx) = real(grid%t_init(i,j,k), kind=c_double)
          end do
        end do
      end do
    end if
    
  end function wrfda_extract_additional_fields

  !> @brief Read and allocate observations using WRFDA's standard pipeline
  !> @details Reads observations from BUFR file and allocates WRFDA iv_type and y_type structures.
  !>          This function calls WRFDA's proven observation reading workflow:
  !>            1. da_setup_obs_structures_bufr → da_read_obs_bufr (reads PREPBUFR)
  !>            2. da_allocate_observations (allocates iv arrays)
  !>            3. da_allocate_y (allocates ob/y arrays)
  !>          CRITICAL: This does NOT compute innovations (iv%synop%u%inv remains empty).
  !>          Innovation computation happens later when da_get_innov_vector is called
  !>          by the observation operator with a background state.
  !> @param[in] grid_ptr C pointer to WRFDA grid structure (for domain bounds)
  !> @param[in] ob_filename C string containing path to BUFR observation file
  !> @param[in] ob_filename_len Length of filename string
  !> @param[out] iv_ptr Output pointer to allocated iv_type structure
  !> @param[out] ob_ptr Output pointer to allocated y_type structure
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_read_and_allocate_observations( &
      grid_ptr, ob_filename, ob_filename_len, iv_ptr, ob_ptr) &
      bind(C, name="wrfda_read_and_allocate_observations")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_char, c_int, c_null_char, c_associated, c_loc
  use module_domain, only: domain, head_grid
  use da_define_structures, only: iv_type, y_type, j_type
  use da_setup_structures, only: da_setup_obs_structures
  use da_wrf_interfaces, only: wrf_message
    implicit none
    
    ! Input parameters
    type(c_ptr), value, intent(in) :: grid_ptr
    type(c_ptr), value, intent(in) :: ob_filename
    integer(c_int), value, intent(in) :: ob_filename_len
    
    ! Output parameters
    type(c_ptr), intent(out) :: iv_ptr
    type(c_ptr), intent(out) :: ob_ptr
    
    type(j_type), save :: wrfda_j
    type(domain), pointer :: grid
    character(len=512) :: bufr_filename
    character(kind=c_char), pointer :: fptr(:)
    integer :: i, max_len
    
    ! Initialize return code
    wrfda_read_and_allocate_observations = 0
    
    ! Validate grid pointer
    if (.not. c_associated(grid_ptr)) then
      call wrf_message("ERROR: grid_ptr is null in wrfda_read_and_allocate_observations")
      wrfda_read_and_allocate_observations = -1
      return
    end if
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(grid_ptr, grid)
    
    ! Convert C filename string to Fortran string
    if (ob_filename_len < 0 .or. ob_filename_len > 1024) then
      call wrf_message("ERROR: Invalid ob_filename_len")
      wrfda_read_and_allocate_observations = -1
      return
    end if
    
    bufr_filename = ""
    call c_f_pointer(ob_filename, fptr, [ob_filename_len])
    max_len = min(ob_filename_len, len(bufr_filename))
    do i = 1, max_len
      if (fptr(i) == c_null_char) exit
      bufr_filename(i:i) = fptr(i)
    end do
    
    
    ! Note: da_setup_obs_structures reads observations based on ob_format in namelist:
    !   ob_format=1 → reads ob.bufr (PREPBUFR format)
    !   ob_format=2 → reads ob.ascii (ASCII format)
    ! File must be in current working directory with these exact names
    
    ! First-time initialization: Set up MPI-like variables
    ! Following da_wrfvar_init1.inc logic for non-DM_PARALLEL builds
    if (.not. wrfda_obs_allocated) then
      num_procs = 1
      myproc = 0
      comm = 0
      rootproc = .true.
      wrfda_obs_allocated = .true.
    end if
    
    ! Call WRFDA's standard observation reading pipeline (the wrapper)
    ! This handles ALL initialization (time_slots, thinning_grid_conv, etc.)
    ! and calls the appropriate reader based on ob_format in namelist
    ! Note: This does NOT compute innovations (grid%xb not needed yet)
    ! With intent(inout), module-level variables work correctly
    call da_setup_obs_structures(head_grid, wrfda_ob, wrfda_iv, wrfda_j)
    
    ! Return C pointers to module-level structures
    iv_ptr = c_loc(wrfda_iv)
    ob_ptr = c_loc(wrfda_ob)
    
  end function wrfda_read_and_allocate_observations

  !> @brief Extract observation counts from WRFDA iv_type structure
  !> @details Returns the number of observations for each observation type
  !> @param[in] iv_ptr C pointer to WRFDA iv_type structure
  !> @param[out] obs_counts Array of observation counts (size: num_ob_indexes)
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_get_obs_type_counts(iv_ptr, obs_counts) &
      bind(C, name="wrfda_get_obs_type_counts")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int
    use da_define_structures, only: iv_type
    use da_control, only: num_ob_indexes
    implicit none
    
    type(c_ptr), value, intent(in) :: iv_ptr
    integer(c_int), intent(out) :: obs_counts(num_ob_indexes)
    
    type(iv_type), pointer :: iv
    integer :: i
    
    wrfda_get_obs_type_counts = 0
    
    if (.not. c_associated(iv_ptr)) then
      obs_counts(:) = 0
      wrfda_get_obs_type_counts = -1
      return
    end if
    
    call c_f_pointer(iv_ptr, iv)
    
    ! Extract observation counts for all types
    ! NOTE: WRFDA observation type indices are named constants (sound=1, synop=2, etc.)
    ! and the info array is indexed by these constants, NOT by 1..num_ob_indexes
    do i = 1, num_ob_indexes
      ! Safety check: verify we're accessing valid memory
      if (i <= size(iv%info)) then
        obs_counts(i) = iv%info(i)%nlocal
      else
        obs_counts(i) = 0
      end if
    end do
    
  end function wrfda_get_obs_type_counts

  !> @brief Extract total observation count from WRFDA iv_type structure
  !> @param[in] iv_ptr C pointer to WRFDA iv_type structure
  !> @param[out] total_count Total number of observations across all types
  !> @return Integer status code (0 = success, non-zero = error)
  integer(c_int) function wrfda_get_total_obs_count(iv_ptr, total_count) &
      bind(C, name="wrfda_get_total_obs_count")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_int
    use da_define_structures, only: iv_type
    use da_control, only: num_ob_indexes
    implicit none
    
    type(c_ptr), value, intent(in) :: iv_ptr
    integer(c_int), intent(out) :: total_count
    
    type(iv_type), pointer :: iv
    integer :: i
    
    wrfda_get_total_obs_count = 0
    total_count = 0
    
    if (.not. c_associated(iv_ptr)) then
      wrfda_get_total_obs_count = -1
      return
    end if
    
    call c_f_pointer(iv_ptr, iv)
    
    ! Sum observation counts across all types
    total_count = 0
    do i = 1, num_ob_indexes
      total_count = total_count + iv%info(i)%nlocal
    end do
    
  end function wrfda_get_total_obs_count

  !> @brief Get pointer to WRFDA iv_type structure (innovation vector)
  !> @details Returns the C pointer to the module-level wrfda_iv structure
  !>          that was allocated during wrfda_read_and_allocate_observations
  !> @return C pointer to iv_type structure, or c_null_ptr if not allocated
  type(c_ptr) function wrfda_get_iv_ptr() bind(C, name="wrfda_get_iv_ptr")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_null_ptr
    use da_define_structures, only: iv_type
    implicit none
    
    if (wrfda_obs_allocated) then
      wrfda_get_iv_ptr = c_loc(wrfda_iv)
    else
      wrfda_get_iv_ptr = c_null_ptr
    end if
  end function wrfda_get_iv_ptr

  !> @brief Get pointer to WRFDA y_type structure (observation values)
  !> @details Returns the C pointer to the module-level wrfda_ob structure
  !>          that was allocated during wrfda_read_and_allocate_observations
  !> @return C pointer to y_type structure, or c_null_ptr if not allocated
  type(c_ptr) function wrfda_get_y_ptr() bind(C, name="wrfda_get_y_ptr")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_null_ptr
    use da_define_structures, only: y_type
    implicit none
    
    if (wrfda_obs_allocated) then
      wrfda_get_y_ptr = c_loc(wrfda_ob)
    else
      wrfda_get_y_ptr = c_null_ptr
    end if
  end function wrfda_get_y_ptr
  
  !> @brief Get pointer to tangent linear y_type (temporary output structure)
  !> @details Returns pointer to wrfda_y_tl which contains H'·δx output
  !>          from da_transform_xtoy, preserving wrfda_ob (observation values)
  !> @return C pointer to wrfda_y_tl, or null if not allocated
  type(c_ptr) function wrfda_get_y_tl_ptr() bind(C, name="wrfda_get_y_tl_ptr")
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_null_ptr
    use da_define_structures, only: y_type
    implicit none
    
    if (wrfda_y_tl_allocated .and. associated(wrfda_y_tl)) then
      wrfda_get_y_tl_ptr = c_loc(wrfda_y_tl)
    else
      wrfda_get_y_tl_ptr = c_null_ptr
    end if
  end function wrfda_get_y_tl_ptr

  !> @brief Transfer WRF fields to background state structure (grid%xb)
  !!
  !! @details This function wraps WRFDA's da_transfer_wrftoxb routine which:
  !! - Transfers WRF native fields (u, v, t, q, etc.) to grid%xb structure
  !! - Computes derived fields (pressure, height, etc.)
  !! - Applies coordinate transformations for Arakawa-C grid
  !! - Prepares grid%xb for use by observation operators
  !!
  !! This MUST be called before using any WRFDA observation operators to ensure
  !! grid%xb is properly populated from the current state.
  !!
  !! @return 0 on success, non-zero on error
  !!
  !! @see da_transfer_wrftoxb.inc in WRFDA source
  function wrfda_transfer_wrftoxb() result(error_code) bind(C, name="wrfda_transfer_wrftoxb")
    use module_domain, only: head_grid, domain
    use module_configure, only: grid_config_rec_type
    use da_define_structures, only: xbx_type
    use da_transfer_model, only: da_transfer_wrftoxb
    implicit none
    
    integer(c_int) :: error_code
    type(domain), pointer :: grid
    type(grid_config_rec_type) :: config_flags
    type(xbx_type) :: xbx
    
    error_code = 0_c_int
    
    ! Get head grid pointer
    grid => head_grid
    if (.not. associated(grid)) then
      call wrf_message("ERROR: head_grid not associated in wrfda_transfer_wrftoxb")
      error_code = -1_c_int
      return
    end if
    
    ! Note: config_flags is not actually used in da_transfer_wrftoxb, 
    ! but is required by the function signature. We initialize it to default values.
    ! (Verified by checking da_transfer_wrftoxb.inc - no config_flags% references)
    
    ! Call WRFDA's da_transfer_wrftoxb to populate grid%xb from WRF fields
    ! This is the standard WRFDA workflow step that must happen before
    ! using observation operators
    call da_transfer_wrftoxb(xbx, grid, config_flags)
    
  end function wrfda_transfer_wrftoxb

  !============================================================================
  ! State Operations - Zero Increment
  !============================================================================
  
  !> @brief Zero the analysis increment (xa) using WRFDA's da_zero_x
  !> @details This function wraps WRFDA's proven da_zero_x subroutine which
  !> zeros all fields in the analysis increment structure. This is used to
  !> initialize the increment at the start of each outer loop iteration.
  !> @param[in] grid_ptr C pointer to domain structure
  !> @return 0 on success, non-zero on error
  integer(c_int) function wrfda_zero_xa(grid_ptr) bind(C, name="wrfda_zero_xa")
    use da_define_structures, only: da_zero_x
    implicit none
    
    type(c_ptr), value :: grid_ptr
    type(domain), pointer :: grid
    
    ! Convert C pointer to Fortran pointer
    call c_f_pointer(grid_ptr, grid)
    
    if (.not. associated(grid)) then
      call wrf_message("ERROR: grid not associated in wrfda_zero_xa")
      wrfda_zero_xa = -1_c_int
      return
    end if
    
    ! Call WRFDA's da_zero_x to zero all increment fields
    ! Note: xa is a member of the grid structure, not a pointer
    call da_zero_x(grid%xa)
    
    wrfda_zero_xa = 0_c_int
    
  end function wrfda_zero_xa

  ! Extract analysis increments from grid%xa to flat arrays
  !> @brief Extract all analysis increment fields to arrays
  !> @details Extracts all 18 identified fields from grid%xa:
  !>          - 3D fields: u, v, w, t, p, q, qt, rh, rho, geoh, wh, qcw, qrn, qci, qsn, qgr
  !>          - 2D fields: psfc, mu
  !> @param[in] grid_ptr Pointer to WRFDA domain structure
  !> @param[out] u, v, w, t, p, q, qt, rh, rho, geoh, wh, qcw, qrn, qci, qsn, qgr 3D arrays (size: nx*ny*nz)
  !> @param[out] psfc, mu 2D arrays (size: nx*ny)
  !> @param[out] nx, ny, nz Grid dimensions
  !> @return 0 on success, -1 on error
  integer(c_int) function wrfda_extract_all_analysis_increments( &
      grid_ptr, &
      u, v, w, t, p, q, qt, rh, rho, geoh, wh, qcw, qrn, qci, qsn, qgr, &
      psfc, mu, &
      nx, ny, nz) &
      bind(C, name="wrfda_extract_all_analysis_increments")
    implicit none
    type(c_ptr), value :: grid_ptr
    real(c_double), intent(out) :: u(*), v(*), w(*), t(*), p(*), q(*)
    real(c_double), intent(out) :: qt(*), rh(*), rho(*), geoh(*), wh(*)
    real(c_double), intent(out) :: qcw(*), qrn(*), qci(*), qsn(*), qgr(*)
    real(c_double), intent(out) :: psfc(*), mu(*)
    integer(c_int), intent(out) :: nx, ny, nz
    
    type(domain), pointer :: grid
    integer :: i, j, k, idx
    
    call c_f_pointer(grid_ptr, grid)
    
    if (.not. associated(grid)) then
      wrfda_extract_all_analysis_increments = -1_c_int
      return
    end if
    
    ! Get dimensions
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    
    ! Extract 3D fields from grid%xa to flat arrays (column-major for C++)
    ! Note: Some fields may be dummy arrays (1x1x1) if not used (e.g., cloud_cv_options <= 1)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          idx = i + (j-1)*nx + (k-1)*nx*ny
          
          ! Primary state variables (always allocated with full dimensions)
          if (associated(grid%xa%u)) u(idx) = real(grid%xa%u(i,j,k), kind=8)
          if (associated(grid%xa%v)) v(idx) = real(grid%xa%v(i,j,k), kind=8)
          if (associated(grid%xa%w)) w(idx) = real(grid%xa%w(i,j,k), kind=8)
          if (associated(grid%xa%t)) t(idx) = real(grid%xa%t(i,j,k), kind=8)
          if (associated(grid%xa%p)) p(idx) = real(grid%xa%p(i,j,k), kind=8)
          if (associated(grid%xa%q)) q(idx) = real(grid%xa%q(i,j,k), kind=8)
          
          ! Additional moisture/thermodynamic variables (check size - may be dummy arrays)
          if (associated(grid%xa%qt)) then
            if (size(grid%xa%qt, 1) >= nx .and. size(grid%xa%qt, 2) >= ny .and. size(grid%xa%qt, 3) >= nz) then
              qt(idx) = real(grid%xa%qt(i,j,k), kind=8)
            end if
          end if
          if (associated(grid%xa%rh)) then
            if (size(grid%xa%rh, 1) >= nx .and. size(grid%xa%rh, 2) >= ny .and. size(grid%xa%rh, 3) >= nz) then
              rh(idx) = real(grid%xa%rh(i,j,k), kind=8)
            end if
          end if
          
          ! Derived fields (check size - may also be conditionally allocated)
          if (associated(grid%xa%rho)) then
            if (size(grid%xa%rho, 1) >= nx .and. size(grid%xa%rho, 2) >= ny .and. size(grid%xa%rho, 3) >= nz) then
              rho(idx) = real(grid%xa%rho(i,j,k), kind=8)
            end if
          end if
          if (associated(grid%xa%geoh)) then
            if (size(grid%xa%geoh, 1) >= nx .and. size(grid%xa%geoh, 2) >= ny .and. size(grid%xa%geoh, 3) >= nz) then
              geoh(idx) = real(grid%xa%geoh(i,j,k), kind=8)
            end if
          end if
          if (associated(grid%xa%wh)) then
            if (size(grid%xa%wh, 1) >= nx .and. size(grid%xa%wh, 2) >= ny .and. size(grid%xa%wh, 3) >= nz) then
              wh(idx) = real(grid%xa%wh(i,j,k), kind=8)
            end if
          end if
          
          ! Hydrometeor fields (check size - conditionally allocated based on cloud_cv_options)
          if (associated(grid%xa%qcw)) then
            if (size(grid%xa%qcw, 1) >= nx .and. size(grid%xa%qcw, 2) >= ny .and. size(grid%xa%qcw, 3) >= nz) then
              qcw(idx) = real(grid%xa%qcw(i,j,k), kind=8)
            end if
          end if
          if (associated(grid%xa%qrn)) then
            if (size(grid%xa%qrn, 1) >= nx .and. size(grid%xa%qrn, 2) >= ny .and. size(grid%xa%qrn, 3) >= nz) then
              qrn(idx) = real(grid%xa%qrn(i,j,k), kind=8)
            end if
          end if
          if (associated(grid%xa%qci)) then
            if (size(grid%xa%qci, 1) >= nx .and. size(grid%xa%qci, 2) >= ny .and. size(grid%xa%qci, 3) >= nz) then
              qci(idx) = real(grid%xa%qci(i,j,k), kind=8)
            end if
          end if
          if (associated(grid%xa%qsn)) then
            if (size(grid%xa%qsn, 1) >= nx .and. size(grid%xa%qsn, 2) >= ny .and. size(grid%xa%qsn, 3) >= nz) then
              qsn(idx) = real(grid%xa%qsn(i,j,k), kind=8)
            end if
          end if
          if (associated(grid%xa%qgr)) then
            if (size(grid%xa%qgr, 1) >= nx .and. size(grid%xa%qgr, 2) >= ny .and. size(grid%xa%qgr, 3) >= nz) then
              qgr(idx) = real(grid%xa%qgr(i,j,k), kind=8)
            end if
          end if
        end do
      end do
    end do
    
    ! Extract 2D surface fields (check size for conditionally allocated fields)
    do j = 1, ny
      do i = 1, nx
        idx = i + (j-1)*nx
        if (associated(grid%xa%psfc)) then
          if (size(grid%xa%psfc, 1) >= nx .and. size(grid%xa%psfc, 2) >= ny) then
            psfc(idx) = real(grid%xa%psfc(i,j), kind=8)
          end if
        end if
        if (associated(grid%xa%mu)) then
          if (size(grid%xa%mu, 1) >= nx .and. size(grid%xa%mu, 2) >= ny) then
            mu(idx) = real(grid%xa%mu(i,j), kind=8)
          end if
        end if
      end do
    end do
    
    wrfda_extract_all_analysis_increments = 0_c_int
    
  end function wrfda_extract_all_analysis_increments

  !> @brief Compute L2 norm of analysis increments
  !> @details Computes ||xa|| = sqrt(sum of squares of all increment fields)
  !>          Includes all meteorological and hydrometeor fields in grid%xa
  !> @param[in] grid_ptr Pointer to WRFDA domain structure
  !> @return L2 norm of analysis increment vector
  function wrfda_xa_norm(grid_ptr) bind(C, name="wrfda_xa_norm") result(norm_val)
    implicit none
    type(c_ptr), value :: grid_ptr
    real(c_double) :: norm_val
    type(domain), pointer :: grid
    integer :: i, j, k
    real(kind=8) :: sum_sq
    integer :: nx, ny, nz
    
    call c_f_pointer(grid_ptr, grid)
    
    if (.not. associated(grid)) then
      norm_val = 0.0_c_double
      return
    end if
    
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    
    sum_sq = 0.0_8
    
    ! Sum squares of all 3D meteorological fields
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! Primary state variables (always allocated with full dimensions)
          if (associated(grid%xa%u)) sum_sq = sum_sq + grid%xa%u(i,j,k)**2
          if (associated(grid%xa%v)) sum_sq = sum_sq + grid%xa%v(i,j,k)**2
          if (associated(grid%xa%w)) sum_sq = sum_sq + grid%xa%w(i,j,k)**2
          if (associated(grid%xa%t)) sum_sq = sum_sq + grid%xa%t(i,j,k)**2
          if (associated(grid%xa%p)) sum_sq = sum_sq + grid%xa%p(i,j,k)**2
          if (associated(grid%xa%q)) sum_sq = sum_sq + grid%xa%q(i,j,k)**2
          
          ! Additional moisture/thermodynamic variables (check size - may be dummy arrays)
          if (associated(grid%xa%qt)) then
            if (size(grid%xa%qt, 1) >= nx .and. size(grid%xa%qt, 2) >= ny .and. size(grid%xa%qt, 3) >= nz) then
              sum_sq = sum_sq + grid%xa%qt(i,j,k)**2
            end if
          end if
          if (associated(grid%xa%rh)) then
            if (size(grid%xa%rh, 1) >= nx .and. size(grid%xa%rh, 2) >= ny .and. size(grid%xa%rh, 3) >= nz) then
              sum_sq = sum_sq + grid%xa%rh(i,j,k)**2
            end if
          end if
          
          ! Derived fields (check size - may also be conditionally allocated)
          if (associated(grid%xa%rho)) then
            if (size(grid%xa%rho, 1) >= nx .and. size(grid%xa%rho, 2) >= ny .and. size(grid%xa%rho, 3) >= nz) then
              sum_sq = sum_sq + grid%xa%rho(i,j,k)**2
            end if
          end if
          if (associated(grid%xa%geoh)) then
            if (size(grid%xa%geoh, 1) >= nx .and. size(grid%xa%geoh, 2) >= ny .and. size(grid%xa%geoh, 3) >= nz) then
              sum_sq = sum_sq + grid%xa%geoh(i,j,k)**2
            end if
          end if
          if (associated(grid%xa%wh)) then
            if (size(grid%xa%wh, 1) >= nx .and. size(grid%xa%wh, 2) >= ny .and. size(grid%xa%wh, 3) >= nz) then
              sum_sq = sum_sq + grid%xa%wh(i,j,k)**2
            end if
          end if
          
          ! Hydrometeor fields (check size - conditionally allocated based on cloud_cv_options)
          if (associated(grid%xa%qcw)) then
            if (size(grid%xa%qcw, 1) >= nx .and. size(grid%xa%qcw, 2) >= ny .and. size(grid%xa%qcw, 3) >= nz) then
              sum_sq = sum_sq + grid%xa%qcw(i,j,k)**2
            end if
          end if
          if (associated(grid%xa%qrn)) then
            if (size(grid%xa%qrn, 1) >= nx .and. size(grid%xa%qrn, 2) >= ny .and. size(grid%xa%qrn, 3) >= nz) then
              sum_sq = sum_sq + grid%xa%qrn(i,j,k)**2
            end if
          end if
          if (associated(grid%xa%qci)) then
            if (size(grid%xa%qci, 1) >= nx .and. size(grid%xa%qci, 2) >= ny .and. size(grid%xa%qci, 3) >= nz) then
              sum_sq = sum_sq + grid%xa%qci(i,j,k)**2
            end if
          end if
          if (associated(grid%xa%qsn)) then
            if (size(grid%xa%qsn, 1) >= nx .and. size(grid%xa%qsn, 2) >= ny .and. size(grid%xa%qsn, 3) >= nz) then
              sum_sq = sum_sq + grid%xa%qsn(i,j,k)**2
            end if
          end if
          if (associated(grid%xa%qgr)) then
            if (size(grid%xa%qgr, 1) >= nx .and. size(grid%xa%qgr, 2) >= ny .and. size(grid%xa%qgr, 3) >= nz) then
              sum_sq = sum_sq + grid%xa%qgr(i,j,k)**2
            end if
          end if
        end do
      end do
    end do
    
    ! Sum squares of 2D surface fields (check size for conditionally allocated fields)
    do j = 1, ny
      do i = 1, nx
        if (associated(grid%xa%psfc)) then
          if (size(grid%xa%psfc, 1) >= nx .and. size(grid%xa%psfc, 2) >= ny) then
            sum_sq = sum_sq + grid%xa%psfc(i,j)**2
          end if
        end if
        if (associated(grid%xa%mu)) then
          if (size(grid%xa%mu, 1) >= nx .and. size(grid%xa%mu, 2) >= ny) then
            sum_sq = sum_sq + grid%xa%mu(i,j)**2
          end if
        end if
      end do
    end do
    
    norm_val = sqrt(sum_sq)
    
  end function wrfda_xa_norm

  !> @brief Compute dot product of two grid%xa structures
  !> @details Computes <xa1, xa2> = sum of element-wise products of all increment fields
  !>          Includes all meteorological and hydrometeor fields in grid%xa
  !> @param[in] grid_ptr1 Pointer to first WRFDA domain structure
  !> @param[in] grid_ptr2 Pointer to second WRFDA domain structure
  !> @return Dot product of two analysis increment vectors
  function wrfda_xa_dot(grid_ptr1, grid_ptr2) bind(C, name="wrfda_xa_dot") result(dot_val)
    implicit none
    type(c_ptr), value :: grid_ptr1, grid_ptr2
    real(c_double) :: dot_val
    type(domain), pointer :: grid1, grid2
    integer :: i, j, k
    real(kind=8) :: sum_prod
    integer :: nx, ny, nz
    
    call c_f_pointer(grid_ptr1, grid1)
    call c_f_pointer(grid_ptr2, grid2)
    
    if (.not. associated(grid1) .or. .not. associated(grid2)) then
      dot_val = 0.0_c_double
      return
    end if
    
    nx = grid1%xp%ide - grid1%xp%ids + 1
    ny = grid1%xp%jde - grid1%xp%jds + 1
    nz = grid1%xp%kde - grid1%xp%kds + 1
    
    sum_prod = 0.0_8
    
    ! Dot product of all 3D meteorological fields
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! Primary state variables (always allocated with full dimensions)
          if (associated(grid1%xa%u) .and. associated(grid2%xa%u)) &
            sum_prod = sum_prod + grid1%xa%u(i,j,k) * grid2%xa%u(i,j,k)
          if (associated(grid1%xa%v) .and. associated(grid2%xa%v)) &
            sum_prod = sum_prod + grid1%xa%v(i,j,k) * grid2%xa%v(i,j,k)
          if (associated(grid1%xa%w) .and. associated(grid2%xa%w)) &
            sum_prod = sum_prod + grid1%xa%w(i,j,k) * grid2%xa%w(i,j,k)
          if (associated(grid1%xa%t) .and. associated(grid2%xa%t)) &
            sum_prod = sum_prod + grid1%xa%t(i,j,k) * grid2%xa%t(i,j,k)
          if (associated(grid1%xa%p) .and. associated(grid2%xa%p)) &
            sum_prod = sum_prod + grid1%xa%p(i,j,k) * grid2%xa%p(i,j,k)
          if (associated(grid1%xa%q) .and. associated(grid2%xa%q)) &
            sum_prod = sum_prod + grid1%xa%q(i,j,k) * grid2%xa%q(i,j,k)
          
          ! Additional moisture/thermodynamic variables (check size - may be dummy arrays)
          if (associated(grid1%xa%qt) .and. associated(grid2%xa%qt)) then
            if (size(grid1%xa%qt, 1) >= nx .and. size(grid1%xa%qt, 2) >= ny .and. size(grid1%xa%qt, 3) >= nz .and. &
                size(grid2%xa%qt, 1) >= nx .and. size(grid2%xa%qt, 2) >= ny .and. size(grid2%xa%qt, 3) >= nz) then
              sum_prod = sum_prod + grid1%xa%qt(i,j,k) * grid2%xa%qt(i,j,k)
            end if
          end if
          if (associated(grid1%xa%rh) .and. associated(grid2%xa%rh)) then
            if (size(grid1%xa%rh, 1) >= nx .and. size(grid1%xa%rh, 2) >= ny .and. size(grid1%xa%rh, 3) >= nz .and. &
                size(grid2%xa%rh, 1) >= nx .and. size(grid2%xa%rh, 2) >= ny .and. size(grid2%xa%rh, 3) >= nz) then
              sum_prod = sum_prod + grid1%xa%rh(i,j,k) * grid2%xa%rh(i,j,k)
            end if
          end if
          
          ! Derived fields (check size - may also be conditionally allocated)
          if (associated(grid1%xa%rho) .and. associated(grid2%xa%rho)) then
            if (size(grid1%xa%rho, 1) >= nx .and. size(grid1%xa%rho, 2) >= ny .and. size(grid1%xa%rho, 3) >= nz .and. &
                size(grid2%xa%rho, 1) >= nx .and. size(grid2%xa%rho, 2) >= ny .and. size(grid2%xa%rho, 3) >= nz) then
              sum_prod = sum_prod + grid1%xa%rho(i,j,k) * grid2%xa%rho(i,j,k)
            end if
          end if
          if (associated(grid1%xa%geoh) .and. associated(grid2%xa%geoh)) then
            if (size(grid1%xa%geoh, 1) >= nx .and. size(grid1%xa%geoh, 2) >= ny .and. size(grid1%xa%geoh, 3) >= nz .and. &
                size(grid2%xa%geoh, 1) >= nx .and. size(grid2%xa%geoh, 2) >= ny .and. size(grid2%xa%geoh, 3) >= nz) then
              sum_prod = sum_prod + grid1%xa%geoh(i,j,k) * grid2%xa%geoh(i,j,k)
            end if
          end if
          if (associated(grid1%xa%wh) .and. associated(grid2%xa%wh)) then
            if (size(grid1%xa%wh, 1) >= nx .and. size(grid1%xa%wh, 2) >= ny .and. size(grid1%xa%wh, 3) >= nz .and. &
                size(grid2%xa%wh, 1) >= nx .and. size(grid2%xa%wh, 2) >= ny .and. size(grid2%xa%wh, 3) >= nz) then
              sum_prod = sum_prod + grid1%xa%wh(i,j,k) * grid2%xa%wh(i,j,k)
            end if
          end if
          
          ! Hydrometeor fields (check size - conditionally allocated based on cloud_cv_options)
          if (associated(grid1%xa%qcw) .and. associated(grid2%xa%qcw)) then
            if (size(grid1%xa%qcw, 1) >= nx .and. size(grid1%xa%qcw, 2) >= ny .and. size(grid1%xa%qcw, 3) >= nz .and. &
                size(grid2%xa%qcw, 1) >= nx .and. size(grid2%xa%qcw, 2) >= ny .and. size(grid2%xa%qcw, 3) >= nz) then
              sum_prod = sum_prod + grid1%xa%qcw(i,j,k) * grid2%xa%qcw(i,j,k)
            end if
          end if
          if (associated(grid1%xa%qrn) .and. associated(grid2%xa%qrn)) then
            if (size(grid1%xa%qrn, 1) >= nx .and. size(grid1%xa%qrn, 2) >= ny .and. size(grid1%xa%qrn, 3) >= nz .and. &
                size(grid2%xa%qrn, 1) >= nx .and. size(grid2%xa%qrn, 2) >= ny .and. size(grid2%xa%qrn, 3) >= nz) then
              sum_prod = sum_prod + grid1%xa%qrn(i,j,k) * grid2%xa%qrn(i,j,k)
            end if
          end if
          if (associated(grid1%xa%qci) .and. associated(grid2%xa%qci)) then
            if (size(grid1%xa%qci, 1) >= nx .and. size(grid1%xa%qci, 2) >= ny .and. size(grid1%xa%qci, 3) >= nz .and. &
                size(grid2%xa%qci, 1) >= nx .and. size(grid2%xa%qci, 2) >= ny .and. size(grid2%xa%qci, 3) >= nz) then
              sum_prod = sum_prod + grid1%xa%qci(i,j,k) * grid2%xa%qci(i,j,k)
            end if
          end if
          if (associated(grid1%xa%qsn) .and. associated(grid2%xa%qsn)) then
            if (size(grid1%xa%qsn, 1) >= nx .and. size(grid1%xa%qsn, 2) >= ny .and. size(grid1%xa%qsn, 3) >= nz .and. &
                size(grid2%xa%qsn, 1) >= nx .and. size(grid2%xa%qsn, 2) >= ny .and. size(grid2%xa%qsn, 3) >= nz) then
              sum_prod = sum_prod + grid1%xa%qsn(i,j,k) * grid2%xa%qsn(i,j,k)
            end if
          end if
          if (associated(grid1%xa%qgr) .and. associated(grid2%xa%qgr)) then
            if (size(grid1%xa%qgr, 1) >= nx .and. size(grid1%xa%qgr, 2) >= ny .and. size(grid1%xa%qgr, 3) >= nz .and. &
                size(grid2%xa%qgr, 1) >= nx .and. size(grid2%xa%qgr, 2) >= ny .and. size(grid2%xa%qgr, 3) >= nz) then
              sum_prod = sum_prod + grid1%xa%qgr(i,j,k) * grid2%xa%qgr(i,j,k)
            end if
          end if
        end do
      end do
    end do
    
    ! Dot product of 2D surface fields (check size for conditionally allocated fields)
    do j = 1, ny
      do i = 1, nx
        if (associated(grid1%xa%psfc) .and. associated(grid2%xa%psfc)) then
          if (size(grid1%xa%psfc, 1) >= nx .and. size(grid1%xa%psfc, 2) >= ny .and. &
              size(grid2%xa%psfc, 1) >= nx .and. size(grid2%xa%psfc, 2) >= ny) then
            sum_prod = sum_prod + grid1%xa%psfc(i,j) * grid2%xa%psfc(i,j)
          end if
        end if
        if (associated(grid1%xa%mu) .and. associated(grid2%xa%mu)) then
          if (size(grid1%xa%mu, 1) >= nx .and. size(grid1%xa%mu, 2) >= ny .and. &
              size(grid2%xa%mu, 1) >= nx .and. size(grid2%xa%mu, 2) >= ny) then
            sum_prod = sum_prod + grid1%xa%mu(i,j) * grid2%xa%mu(i,j)
          end if
        end if
      end do
    end do
    
    dot_val = sum_prod
    
  end function wrfda_xa_dot

  !============================================================================
  ! Increment Vector Operations
  !============================================================================
  
  ! AXPY operation on analysis increments: grid1%xa = grid1%xa + alpha * grid2%xa
  !> @brief Perform AXPY operation on grid%xa structures
  !> @details Computes xa1 = xa1 + alpha * xa2 for all increment fields
  !>          Includes all meteorological and hydrometeor fields in grid%xa
  !> @param[in] alpha Scalar multiplier
  !> @param[in,out] grid_ptr1 Pointer to first WRFDA domain (updated in-place)
  !> @param[in] grid_ptr2 Pointer to second WRFDA domain (read-only)
  subroutine wrfda_xa_axpy(alpha, grid_ptr1, grid_ptr2) bind(C, name="wrfda_xa_axpy")
    implicit none
    real(c_double), value :: alpha
    type(c_ptr), value :: grid_ptr1, grid_ptr2
    type(domain), pointer :: grid1, grid2
    integer :: i, j, k, nx, ny, nz
    real(kind=4) :: alpha_real
    
    call c_f_pointer(grid_ptr1, grid1)
    call c_f_pointer(grid_ptr2, grid2)
    
    if (.not. associated(grid1) .or. .not. associated(grid2)) return
    
    nx = grid1%xp%ide - grid1%xp%ids + 1
    ny = grid1%xp%jde - grid1%xp%jds + 1
    nz = grid1%xp%kde - grid1%xp%kds + 1
    alpha_real = real(alpha, kind=4)
    
    ! AXPY for all 3D meteorological fields
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! Primary state variables (always allocated with full dimensions)
          if (associated(grid1%xa%u) .and. associated(grid2%xa%u)) &
            grid1%xa%u(i,j,k) = grid1%xa%u(i,j,k) + alpha_real * grid2%xa%u(i,j,k)
          if (associated(grid1%xa%v) .and. associated(grid2%xa%v)) &
            grid1%xa%v(i,j,k) = grid1%xa%v(i,j,k) + alpha_real * grid2%xa%v(i,j,k)
          if (associated(grid1%xa%w) .and. associated(grid2%xa%w)) &
            grid1%xa%w(i,j,k) = grid1%xa%w(i,j,k) + alpha_real * grid2%xa%w(i,j,k)
          if (associated(grid1%xa%t) .and. associated(grid2%xa%t)) &
            grid1%xa%t(i,j,k) = grid1%xa%t(i,j,k) + alpha_real * grid2%xa%t(i,j,k)
          if (associated(grid1%xa%p) .and. associated(grid2%xa%p)) &
            grid1%xa%p(i,j,k) = grid1%xa%p(i,j,k) + alpha_real * grid2%xa%p(i,j,k)
          if (associated(grid1%xa%q) .and. associated(grid2%xa%q)) &
            grid1%xa%q(i,j,k) = grid1%xa%q(i,j,k) + alpha_real * grid2%xa%q(i,j,k)
          
          ! Additional moisture/thermodynamic variables (check size - may be dummy arrays)
          if (associated(grid1%xa%qt) .and. associated(grid2%xa%qt)) then
            if (size(grid1%xa%qt, 1) >= nx .and. size(grid1%xa%qt, 2) >= ny .and. size(grid1%xa%qt, 3) >= nz .and. &
                size(grid2%xa%qt, 1) >= nx .and. size(grid2%xa%qt, 2) >= ny .and. size(grid2%xa%qt, 3) >= nz) then
              grid1%xa%qt(i,j,k) = grid1%xa%qt(i,j,k) + alpha_real * grid2%xa%qt(i,j,k)
            end if
          end if
          if (associated(grid1%xa%rh) .and. associated(grid2%xa%rh)) then
            if (size(grid1%xa%rh, 1) >= nx .and. size(grid1%xa%rh, 2) >= ny .and. size(grid1%xa%rh, 3) >= nz .and. &
                size(grid2%xa%rh, 1) >= nx .and. size(grid2%xa%rh, 2) >= ny .and. size(grid2%xa%rh, 3) >= nz) then
              grid1%xa%rh(i,j,k) = grid1%xa%rh(i,j,k) + alpha_real * grid2%xa%rh(i,j,k)
            end if
          end if
          
          ! Derived fields (check size - may also be conditionally allocated)
          if (associated(grid1%xa%rho) .and. associated(grid2%xa%rho)) then
            if (size(grid1%xa%rho, 1) >= nx .and. size(grid1%xa%rho, 2) >= ny .and. size(grid1%xa%rho, 3) >= nz .and. &
                size(grid2%xa%rho, 1) >= nx .and. size(grid2%xa%rho, 2) >= ny .and. size(grid2%xa%rho, 3) >= nz) then
              grid1%xa%rho(i,j,k) = grid1%xa%rho(i,j,k) + alpha_real * grid2%xa%rho(i,j,k)
            end if
          end if
          if (associated(grid1%xa%geoh) .and. associated(grid2%xa%geoh)) then
            if (size(grid1%xa%geoh, 1) >= nx .and. size(grid1%xa%geoh, 2) >= ny .and. size(grid1%xa%geoh, 3) >= nz .and. &
                size(grid2%xa%geoh, 1) >= nx .and. size(grid2%xa%geoh, 2) >= ny .and. size(grid2%xa%geoh, 3) >= nz) then
              grid1%xa%geoh(i,j,k) = grid1%xa%geoh(i,j,k) + alpha_real * grid2%xa%geoh(i,j,k)
            end if
          end if
          if (associated(grid1%xa%wh) .and. associated(grid2%xa%wh)) then
            if (size(grid1%xa%wh, 1) >= nx .and. size(grid1%xa%wh, 2) >= ny .and. size(grid1%xa%wh, 3) >= nz .and. &
                size(grid2%xa%wh, 1) >= nx .and. size(grid2%xa%wh, 2) >= ny .and. size(grid2%xa%wh, 3) >= nz) then
              grid1%xa%wh(i,j,k) = grid1%xa%wh(i,j,k) + alpha_real * grid2%xa%wh(i,j,k)
            end if
          end if
          
          ! Hydrometeor fields (check size - conditionally allocated based on cloud_cv_options)
          if (associated(grid1%xa%qcw) .and. associated(grid2%xa%qcw)) then
            if (size(grid1%xa%qcw, 1) >= nx .and. size(grid1%xa%qcw, 2) >= ny .and. size(grid1%xa%qcw, 3) >= nz .and. &
                size(grid2%xa%qcw, 1) >= nx .and. size(grid2%xa%qcw, 2) >= ny .and. size(grid2%xa%qcw, 3) >= nz) then
              grid1%xa%qcw(i,j,k) = grid1%xa%qcw(i,j,k) + alpha_real * grid2%xa%qcw(i,j,k)
            end if
          end if
          if (associated(grid1%xa%qrn) .and. associated(grid2%xa%qrn)) then
            if (size(grid1%xa%qrn, 1) >= nx .and. size(grid1%xa%qrn, 2) >= ny .and. size(grid1%xa%qrn, 3) >= nz .and. &
                size(grid2%xa%qrn, 1) >= nx .and. size(grid2%xa%qrn, 2) >= ny .and. size(grid2%xa%qrn, 3) >= nz) then
              grid1%xa%qrn(i,j,k) = grid1%xa%qrn(i,j,k) + alpha_real * grid2%xa%qrn(i,j,k)
            end if
          end if
          if (associated(grid1%xa%qci) .and. associated(grid2%xa%qci)) then
            if (size(grid1%xa%qci, 1) >= nx .and. size(grid1%xa%qci, 2) >= ny .and. size(grid1%xa%qci, 3) >= nz .and. &
                size(grid2%xa%qci, 1) >= nx .and. size(grid2%xa%qci, 2) >= ny .and. size(grid2%xa%qci, 3) >= nz) then
              grid1%xa%qci(i,j,k) = grid1%xa%qci(i,j,k) + alpha_real * grid2%xa%qci(i,j,k)
            end if
          end if
          if (associated(grid1%xa%qsn) .and. associated(grid2%xa%qsn)) then
            if (size(grid1%xa%qsn, 1) >= nx .and. size(grid1%xa%qsn, 2) >= ny .and. size(grid1%xa%qsn, 3) >= nz .and. &
                size(grid2%xa%qsn, 1) >= nx .and. size(grid2%xa%qsn, 2) >= ny .and. size(grid2%xa%qsn, 3) >= nz) then
              grid1%xa%qsn(i,j,k) = grid1%xa%qsn(i,j,k) + alpha_real * grid2%xa%qsn(i,j,k)
            end if
          end if
          if (associated(grid1%xa%qgr) .and. associated(grid2%xa%qgr)) then
            if (size(grid1%xa%qgr, 1) >= nx .and. size(grid1%xa%qgr, 2) >= ny .and. size(grid1%xa%qgr, 3) >= nz .and. &
                size(grid2%xa%qgr, 1) >= nx .and. size(grid2%xa%qgr, 2) >= ny .and. size(grid2%xa%qgr, 3) >= nz) then
              grid1%xa%qgr(i,j,k) = grid1%xa%qgr(i,j,k) + alpha_real * grid2%xa%qgr(i,j,k)
            end if
          end if
        end do
      end do
    end do
    
    ! AXPY for 2D surface fields (check size for conditionally allocated fields)
    do j = 1, ny
      do i = 1, nx
        if (associated(grid1%xa%psfc) .and. associated(grid2%xa%psfc)) then
          if (size(grid1%xa%psfc, 1) >= nx .and. size(grid1%xa%psfc, 2) >= ny .and. &
              size(grid2%xa%psfc, 1) >= nx .and. size(grid2%xa%psfc, 2) >= ny) then
            grid1%xa%psfc(i,j) = grid1%xa%psfc(i,j) + alpha_real * grid2%xa%psfc(i,j)
          end if
        end if
        if (associated(grid1%xa%mu) .and. associated(grid2%xa%mu)) then
          if (size(grid1%xa%mu, 1) >= nx .and. size(grid1%xa%mu, 2) >= ny .and. &
              size(grid2%xa%mu, 1) >= nx .and. size(grid2%xa%mu, 2) >= ny) then
            grid1%xa%mu(i,j) = grid1%xa%mu(i,j) + alpha_real * grid2%xa%mu(i,j)
          end if
        end if
      end do
    end do
    
  end subroutine wrfda_xa_axpy

  !> @brief Scale analysis increments by a scalar
  !> @details Computes xa = alpha * xa for all increment fields
  !>          Includes all meteorological and hydrometeor fields in grid%xa
  !> @param[in] alpha Scalar multiplier
  !> @param[in,out] grid_ptr Pointer to WRFDA domain (updated in-place)
  subroutine wrfda_xa_scale(alpha, grid_ptr) bind(C, name="wrfda_xa_scale")
    implicit none
    real(c_double), value :: alpha
    type(c_ptr), value :: grid_ptr
    type(domain), pointer :: grid
    integer :: i, j, k, nx, ny, nz
    real(kind=4) :: alpha_real
    
    call c_f_pointer(grid_ptr, grid)
    
    if (.not. associated(grid)) return
    
    nx = grid%xp%ide - grid%xp%ids + 1
    ny = grid%xp%jde - grid%xp%jds + 1
    nz = grid%xp%kde - grid%xp%kds + 1
    alpha_real = real(alpha, kind=4)
    
    ! Scale all 3D meteorological fields
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! Primary state variables (always allocated with full dimensions)
          if (associated(grid%xa%u)) grid%xa%u(i,j,k) = alpha_real * grid%xa%u(i,j,k)
          if (associated(grid%xa%v)) grid%xa%v(i,j,k) = alpha_real * grid%xa%v(i,j,k)
          if (associated(grid%xa%w)) grid%xa%w(i,j,k) = alpha_real * grid%xa%w(i,j,k)
          if (associated(grid%xa%t)) grid%xa%t(i,j,k) = alpha_real * grid%xa%t(i,j,k)
          if (associated(grid%xa%p)) grid%xa%p(i,j,k) = alpha_real * grid%xa%p(i,j,k)
          if (associated(grid%xa%q)) grid%xa%q(i,j,k) = alpha_real * grid%xa%q(i,j,k)
          
          ! Additional moisture/thermodynamic variables (check size - may be dummy arrays)
          if (associated(grid%xa%qt)) then
            if (size(grid%xa%qt, 1) >= nx .and. size(grid%xa%qt, 2) >= ny .and. size(grid%xa%qt, 3) >= nz) then
              grid%xa%qt(i,j,k) = alpha_real * grid%xa%qt(i,j,k)
            end if
          end if
          if (associated(grid%xa%rh)) then
            if (size(grid%xa%rh, 1) >= nx .and. size(grid%xa%rh, 2) >= ny .and. size(grid%xa%rh, 3) >= nz) then
              grid%xa%rh(i,j,k) = alpha_real * grid%xa%rh(i,j,k)
            end if
          end if
          
          ! Derived fields (check size - may also be conditionally allocated)
          if (associated(grid%xa%rho)) then
            if (size(grid%xa%rho, 1) >= nx .and. size(grid%xa%rho, 2) >= ny .and. size(grid%xa%rho, 3) >= nz) then
              grid%xa%rho(i,j,k) = alpha_real * grid%xa%rho(i,j,k)
            end if
          end if
          if (associated(grid%xa%geoh)) then
            if (size(grid%xa%geoh, 1) >= nx .and. size(grid%xa%geoh, 2) >= ny .and. size(grid%xa%geoh, 3) >= nz) then
              grid%xa%geoh(i,j,k) = alpha_real * grid%xa%geoh(i,j,k)
            end if
          end if
          if (associated(grid%xa%wh)) then
            if (size(grid%xa%wh, 1) >= nx .and. size(grid%xa%wh, 2) >= ny .and. size(grid%xa%wh, 3) >= nz) then
              grid%xa%wh(i,j,k) = alpha_real * grid%xa%wh(i,j,k)
            end if
          end if
          
          ! Hydrometeor fields (check size - conditionally allocated based on cloud_cv_options)
          if (associated(grid%xa%qcw)) then
            if (size(grid%xa%qcw, 1) >= nx .and. size(grid%xa%qcw, 2) >= ny .and. size(grid%xa%qcw, 3) >= nz) then
              grid%xa%qcw(i,j,k) = alpha_real * grid%xa%qcw(i,j,k)
            end if
          end if
          if (associated(grid%xa%qrn)) then
            if (size(grid%xa%qrn, 1) >= nx .and. size(grid%xa%qrn, 2) >= ny .and. size(grid%xa%qrn, 3) >= nz) then
              grid%xa%qrn(i,j,k) = alpha_real * grid%xa%qrn(i,j,k)
            end if
          end if
          if (associated(grid%xa%qci)) then
            if (size(grid%xa%qci, 1) >= nx .and. size(grid%xa%qci, 2) >= ny .and. size(grid%xa%qci, 3) >= nz) then
              grid%xa%qci(i,j,k) = alpha_real * grid%xa%qci(i,j,k)
            end if
          end if
          if (associated(grid%xa%qsn)) then
            if (size(grid%xa%qsn, 1) >= nx .and. size(grid%xa%qsn, 2) >= ny .and. size(grid%xa%qsn, 3) >= nz) then
              grid%xa%qsn(i,j,k) = alpha_real * grid%xa%qsn(i,j,k)
            end if
          end if
          if (associated(grid%xa%qgr)) then
            if (size(grid%xa%qgr, 1) >= nx .and. size(grid%xa%qgr, 2) >= ny .and. size(grid%xa%qgr, 3) >= nz) then
              grid%xa%qgr(i,j,k) = alpha_real * grid%xa%qgr(i,j,k)
            end if
          end if
        end do
      end do
    end do
    
    ! Scale 2D surface fields (check size for conditionally allocated fields)
    do j = 1, ny
      do i = 1, nx
        if (associated(grid%xa%psfc)) then
          if (size(grid%xa%psfc, 1) >= nx .and. size(grid%xa%psfc, 2) >= ny) then
            grid%xa%psfc(i,j) = alpha_real * grid%xa%psfc(i,j)
          end if
        end if
        if (associated(grid%xa%mu)) then
          if (size(grid%xa%mu, 1) >= nx .and. size(grid%xa%mu, 2) >= ny) then
            grid%xa%mu(i,j) = alpha_real * grid%xa%mu(i,j)
          end if
        end if
      end do
    end do
    
  end subroutine wrfda_xa_scale

end module metada_wrfda_dispatch