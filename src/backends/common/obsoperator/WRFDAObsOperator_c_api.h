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

// New function that properly separates total state into background and
// increments
int wrfda_xtoy_apply_grid_with_background(
    const char* operator_family, int nx, int ny, int nz, const double* u_total,
    const double* v_total, const double* t_total, const double* q_total,
    const double* psfc_total,  // total state (xb + δx)
    const double* u_bg, const double* v_bg, const double* t_bg,
    const double* q_bg, const double* psfc_bg,   // background state (xb)
    const double* lats2d, const double* lons2d,  // size nx*ny
    const double* levels,                        // size nz
    int num_obs, const double* obs_lats, const double* obs_lons,
    const double* obs_levels, double* out_y);

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

// Enhanced function for constructing iv_type with detailed observation
// information
int wrfda_construct_iv_from_observations(
    const char* operator_family, int nx, int ny, int nz, const double* lats2d,
    const double* lons2d, const double* levels,  // size nx*ny, nz
    int num_obs, const double* obs_lats, const double* obs_lons,
    const double* obs_levels, const double* obs_values,
    const double* obs_errors,  // size num_obs
    const char* obs_types,
    const char* obs_station_ids,  // size num_obs * string_length
    const double* obs_elevations, const double* obs_pressures,
    const double* obs_heights,                            // size num_obs
    const int* obs_qc_flags, const int* obs_usage_flags,  // size num_obs
    const double* obs_time_offsets,                       // size num_obs
    void* iv_ptr  // pointer to iv_type structure
);

// Main WRFDA innovation vector computation (calls da_get_innov_vector)
// This is the comprehensive driver routine that handles all observation types
int wrfda_get_innov_vector(const int* it, const void* domain_ptr,
                           const void* ob_ptr, const void* iv_ptr,
                           const void* config_flags_ptr);

// Helper function to construct WRFDA domain structure from flat arrays
int wrfda_construct_domain_from_arrays(
    const int* nx, const int* ny, const int* nz, const double* u,
    const double* v, const double* t, const double* q, const double* psfc,
    const double* ph, const double* phb, const double* hf, const double* hgt,
    const double* lats2d, const double* lons2d, const double* levels,
    void** domain_ptr);

// Construct y_type from observation data
void* wrfda_construct_y_type(int* num_obs, int* num_levels,
                             const double* u_values, const double* v_values,
                             const double* t_values, const double* p_values,
                             const double* q_values, const double* u_errors,
                             const double* v_errors, const double* t_errors,
                             const double* p_errors, const double* q_errors,
                             const int* u_available, const int* v_available,
                             const int* t_available, const int* p_available,
                             const int* q_available, const double* lats,
                             const double* lons, const double* levels,
                             const char* obs_types, const char* family);

// Construct iv_type from observation data
void* wrfda_construct_iv_type(
    int* num_obs, int* num_levels, const double* u_values,
    const double* v_values, const double* t_values, const double* p_values,
    const double* q_values, const double* u_errors, const double* v_errors,
    const double* t_errors, const double* p_errors, const double* q_errors,
    const int* u_qc, const int* v_qc, const int* t_qc, const int* p_qc,
    const int* q_qc, const int* u_available, const int* v_available,
    const int* t_available, const int* p_available, const int* q_available,
    const double* lats, const double* lons, const double* levels,
    const double* elevations, const char* obs_types, const char* family,
    void* domain_ptr);

// Construct config_flags for WRFDA
void* wrfda_construct_config_flags();

// Extract innovation values from iv_type structure
int wrfda_extract_innovations(void* iv_ptr, const char* family,
                              double* innovations, int* num_innovations,
                              int* max_innovations);

// Initialize WRFDA variables for 3D-Var analysis
void initialize_wrfda_3dvar();

// Initialize map projection with grid parameters
void initialize_map_projection_c(const int* map_proj, const double* cen_lat,
                                 const double* cen_lon, const double* dx,
                                 const double* stand_lon,
                                 const double* truelat1,
                                 const double* truelat2);

// Initialize WRFDA module-level variables (kts, kte, sfc_assi_options, etc.)
int initialize_wrfda_module_variables(void* domain_ptr);

#ifdef __cplusplus
}
#endif
