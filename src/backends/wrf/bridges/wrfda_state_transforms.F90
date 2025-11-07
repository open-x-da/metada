!> @file wrfda_state_transforms.F90
!> @brief Fortran-C bridges for WRFDA state transformation routines
!> @details Provides C-callable interfaces to WRFDA's proven state transformation
!>          routines without modifying WRFDA source code

module wrfda_state_transforms_bridge
  use iso_c_binding
  use module_domain, only: domain
  use module_configure, only: grid_config_rec_type
  use da_vtox_transforms, only: da_transform_xtoxa
  use da_transfer_model, only: da_transfer_xatowrf
  use da_wrfvar_io, only: da_update_firstguess
  implicit none

  private
  public :: wrfda_transform_xtoxa_c
  public :: wrfda_transfer_xatowrf_c
  public :: wrfda_update_firstguess_c

contains

  !> @brief C-callable bridge to WRFDA's da_transform_xtoxa
  !> @details Computes diagnostic increments (p, rho, geoh) from prognostic increments.
  !>          This routine is called during cost function evaluation to prepare
  !>          diagnostic variables needed by observation operators.
  !>
  !> @param[in] grid_ptr C pointer to WRF grid (domain) structure
  !> @param[out] error_code 0 = success, non-zero = error
  !>
  !> @note This routine should be called after da_transform_vtox and before
  !>       observation operators are applied
  subroutine wrfda_transform_xtoxa_c(grid_ptr, error_code) bind(C, name="wrfda_transform_xtoxa")
    implicit none

    ! Arguments
    type(c_ptr), value, intent(in) :: grid_ptr
    integer(c_int), intent(out) :: error_code

    ! Local variables
    type(domain), pointer :: grid

    ! Initialize error code
    error_code = 0

    ! Safety check for null pointer
    if (.not. c_associated(grid_ptr)) then
      error_code = 1
      return
    end if

    ! Convert C pointer to Fortran pointer
    call c_f_pointer(grid_ptr, grid)

    ! Call WRFDA's proven routine
    ! This computes:
    !   - grid%xa%p (pressure increments)
    !   - grid%xa%rho (density increments)
    !   - grid%xa%geoh (geopotential height increments)
    !   - Derived variables for observation operators
    call da_transform_xtoxa(grid)

  end subroutine wrfda_transform_xtoxa_c

  !> @brief C-callable bridge to WRFDA's da_transfer_xatowrf
  !> @details Adds analysis increments to background state and updates all WRF
  !>          prognostic variables. This routine handles:
  !>          - Conversion of specific humidity to mixing ratio
  !>          - Computation of dry air mass increments
  !>          - Conversion of temperature to potential temperature
  !>          - Computation of geopotential height for WRF's vertical coordinate
  !>          - Arakawa-C grid staggering (A-grid to C-grid conversion)
  !>          - Positivity constraints (e.g., moisture >= 0)
  !>          - Recomputation of 2m/10m diagnostic fields
  !>
  !> @param[in] grid_ptr C pointer to WRF grid (domain) structure
  !> @param[in] config_flags_ptr C pointer to WRF config_flags structure
  !> @param[out] error_code 0 = success, non-zero = error
  !>
  !> @note After this routine, the WRF state contains the full analysis:
  !>       x_analysis = x_background + increment
  !>
  !> @note This routine modifies the following WRF fields:
  !>       - grid%u_2, grid%v_2, grid%w_2 (wind components)
  !>       - grid%mu_2 (dry air mass)
  !>       - grid%ph_2 (geopotential)
  !>       - grid%t_2 (potential temperature perturbation)
  !>       - grid%moist (moisture species)
  !>       - grid%p (perturbation pressure)
  !>       - grid%psfc (surface pressure)
  !>       - grid%t2, grid%q2, grid%u10, grid%v10, grid%th2 (diagnostic fields)
  subroutine wrfda_transfer_xatowrf_c(grid_ptr, config_flags_ptr, error_code) &
      bind(C, name="wrfda_transfer_xatowrf")
    implicit none

    ! Arguments
    type(c_ptr), value, intent(in) :: grid_ptr
    type(c_ptr), value, intent(in) :: config_flags_ptr
    integer(c_int), intent(out) :: error_code

    ! Local variables
    type(domain), pointer :: grid
    type(grid_config_rec_type), pointer :: config_flags

    ! Initialize error code
    error_code = 0

    ! Safety checks for null pointers
    if (.not. c_associated(grid_ptr)) then
      error_code = 1
      return
    end if

    if (.not. c_associated(config_flags_ptr)) then
      error_code = 2
      return
    end if

    ! Convert C pointers to Fortran pointers
    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(config_flags_ptr, config_flags)

    ! Call WRFDA's proven routine
    ! This performs the complete increment-to-state transformation
    ! and updates all WRF prognostic variables
    call da_transfer_xatowrf(grid, config_flags)

  end subroutine wrfda_transfer_xatowrf_c

  !> @brief C-callable bridge to WRFDA's da_update_firstguess
  !> @details Uses WRFDA's battle-tested NetCDF update routine that copies the
  !>          background (fg) file to the specified output file and updates only
  !>          the analysis variables touched during assimilation.
  !>
  !> @param[in] grid_ptr Pointer to WRF grid (domain) structure
  !> @param[in] filename Optional output filename (pass length 0 for default)
  !> @param[in] filename_len Length of filename buffer
  !> @param[out] error_code 0 = success, non-zero = error
  subroutine wrfda_update_firstguess_c(grid_ptr, filename, filename_len, error_code) &
      bind(C, name="wrfda_update_firstguess")
    implicit none

    type(c_ptr), value, intent(in) :: grid_ptr
    character(kind=c_char), dimension(*), intent(in) :: filename
    integer(c_int), value, intent(in) :: filename_len
    integer(c_int), intent(out) :: error_code

    type(domain), pointer :: grid
    character(len=:), allocatable :: fname
    character(len=1) :: ch
    logical :: has_filename
    integer :: i
    integer :: effective_len

    error_code = 0

    if (.not. c_associated(grid_ptr)) then
      error_code = 1
      return
    end if

    if (filename_len < 0) then
      error_code = 3
      return
    end if

    call c_f_pointer(grid_ptr, grid)

    has_filename = .false.
    if (filename_len > 0) then
      effective_len = 0
      do i = 1, filename_len
        if (filename(i) == c_null_char) exit
        effective_len = effective_len + 1
      end do
      if (effective_len > 0) then
        allocate(character(len=effective_len) :: fname)
        do i = 1, effective_len
          ch = transfer(filename(i), ch)
          fname(i:i) = ch
        end do
        fname = trim(adjustl(fname))
        if (len_trim(fname) > 0) then
          has_filename = .true.
        else
          deallocate(fname)
        end if
      end if
    end if

    if (has_filename) then
      call da_update_firstguess(grid, fname)
    else
      call da_update_firstguess(grid)
    end if

    if (allocated(fname)) deallocate(fname)
  end subroutine wrfda_update_firstguess_c

end module wrfda_state_transforms_bridge

