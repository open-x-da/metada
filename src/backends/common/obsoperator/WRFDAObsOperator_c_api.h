/**
 * @file WRFDAObsOperator_c_api.h
 * @brief C-facing API stub for bridging to Fortran WRFDA observation operators
 *
 * @details
 * This header declares a minimal C API intended to be implemented by a small
 * Fortran/C interoperability library that calls WRFDA routines under
 * D:/linux/WRF/var/da (da_transform_xtoy_* and adjoints). The implementation
 * can be provided separately per platform/toolchain and linked externally.
 *
 * The C++ WRFDAObsOperator can optionally dlopen/LoadLibrary such bridge later
 * or link to it directly when available.
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Forward operator: y = H(x)
// Returns 0 on success, non-zero on failure.
int wrfda_xtoy_apply(
    const char* operator_family,        // e.g., "metar", "gpspw"
    const char* wrfda_root,             // path to WRF/var/da
    const double* state_values,         // flattened x array (backend-defined layout)
    int nx, int ny, int nz,             // grid dims
    const double* obs_lats,             // observation lats
    const double* obs_lons,             // observation lons
    const double* obs_levels,           // observation vertical coordinate (pressure/level)
    int num_obs,                        // number of observations
    double* out_y                       // output y array of size num_obs
);

// Adjoint operator: x += H^T(delta_y)
int wrfda_xtoy_adjoint(
    const char* operator_family,
    const char* wrfda_root,
    const double* delta_y,              // size: num_obs
    const double* obs_lats,
    const double* obs_lons,
    const double* obs_levels,
    int num_obs,
    double* inout_state_values,         // flattened x array to accumulate into
    int nx, int ny, int nz
);

// Handle-based API (preferred when WRFDA domain/iv/y are already allocated)
// Opaque pointers are expected to point to WRFDA types:
//   domain :: WRFDA domain/grid structure
//   iv     :: iv_type containing per-family interpolation/index info
//   y/x    :: y_type or x_type for outputs (see function)
// Returns 0 on success
int wrfda_xtoy_apply_handles(
    const char* operator_family,
    const void* domain_ptr,
    const void* iv_ptr,
    void* y_ptr
);

int wrfda_xtoy_adjoint_handles(
    const char* operator_family,
    const void* domain_ptr,
    const void* iv_ptr,
    const void* jo_grad_y_ptr,
    void* jo_grad_x_ptr
);

#ifdef __cplusplus
}
#endif


