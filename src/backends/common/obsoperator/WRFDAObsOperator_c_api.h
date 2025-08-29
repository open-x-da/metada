/**
 * @file WRFDAObsOperator_c_api.h
 * @brief C-facing API for bridging to Fortran WRFDA observation operators
 *
 * @details
 * This header declares a minimal C API implemented by a Fortran/C
 * interoperability library that calls WRFDA routines under D:/linux/WRF/var/da
 * (da_transform_xtoy_* and adjoints). The implementation is provided by the
 * wrfda_bridge library.
 *
 * The C++ WRFDAObsOperator links to this bridge library directly.
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Forward operator: y = H(x)
// Returns 0 on success, non-zero on failure.
int wrfda_xtoy_apply(
    const char* operator_family,  // e.g., "metar", "gpspw"
    const char* wrfda_root,       // path to WRF/var/da
    const double* state_values,   // flattened x array (backend-defined layout)
    int nx, int ny, int nz,       // grid dims
    const double* obs_lats,       // observation lats
    const double* obs_lons,       // observation lons
    const double*
        obs_levels,  // observation vertical coordinate (pressure/level)
    int num_obs,     // number of observations
    double* out_y    // output y array of size num_obs
);

// Adjoint operator: x += H^T(delta_y)
int wrfda_xtoy_adjoint(
    const char* operator_family, const char* wrfda_root,
    const double* delta_y,  // size: num_obs
    const double* obs_lats, const double* obs_lons, const double* obs_levels,
    int num_obs,
    double* inout_state_values,  // flattened x array to accumulate into
    int nx, int ny, int nz);

// Handle-based API (preferred when WRFDA domain/iv/y are already allocated)
// Opaque pointers are expected to point to WRFDA types:
//   domain :: WRFDA domain/grid structure
//   iv     :: iv_type containing per-family interpolation/index info
//   y/x    :: y_type or x_type for outputs (see function)
// Returns 0 on success
int wrfda_xtoy_apply_handles(const char* operator_family,
                             const void* domain_ptr, const void* iv_ptr,
                             void* y_ptr);

int wrfda_xtoy_adjoint_handles(const char* operator_family,
                               const void* domain_ptr, const void* iv_ptr,
                               const void* jo_grad_y_ptr, void* jo_grad_x_ptr);

// Array-based (grid) API for real WRFDA call with per-variable fields and grid
// metadata
int wrfda_xtoy_apply_grid(const char* operator_family, int nx, int ny, int nz,
                          const double* u, const double* v, const double* t,
                          const double* q,
                          const double* psfc,  // size nx*ny
                          const double* lats2d,
                          const double* lons2d,  // size nx*ny
                          const double* levels,  // size nz
                          int num_obs, const double* obs_lats,
                          const double* obs_lons, const double* obs_levels,
                          double* out_y);

int wrfda_xtoy_adjoint_grid(const char* operator_family, int nx, int ny, int nz,
                            const double* delta_y,  // size num_obs
                            const double* lats2d, const double* lons2d,
                            const double* levels,  // size nz
                            int num_obs, const double* obs_lats,
                            const double* obs_lons, const double* obs_levels,
                            double* inout_u, double* inout_v, double* inout_t,
                            double* inout_q,
                            double* inout_psfc  // size nx*ny
);

// Profile-capable array APIs (multi-level per obs, flattened per-obs levels)
int wrfda_xtoy_apply_profiles(
    const char* operator_family,  // e.g., "airep:u", "pilot:u", "sound:t"
    int nx, int ny, int nz, const double* u, const double* v, const double* t,
    const double* q, const double* psfc, const double* lats2d,
    const double* lons2d, const double* levels, int num_obs,
    const int* obs_counts,          // length num_obs; levels per obs
    const double* obs_lats,         // length num_obs
    const double* obs_lons,         // length num_obs
    const double* obs_levels_flat,  // length sum(obs_counts)
    double* out_y_flat              // length sum(obs_counts)
);

int wrfda_xtoy_adjoint_profiles(
    const char* operator_family, int nx, int ny, int nz,
    const double* delta_y_flat,  // length sum(obs_counts)
    const double* lats2d, const double* lons2d, const double* levels,
    int num_obs, const int* obs_counts, const double* obs_lats,
    const double* obs_lons, const double* obs_levels_flat, double* inout_u,
    double* inout_v, double* inout_t, double* inout_q, double* inout_psfc);

#ifdef __cplusplus
}
#endif
