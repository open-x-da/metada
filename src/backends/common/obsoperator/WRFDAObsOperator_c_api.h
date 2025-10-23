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

// Array-based (grid) API for real WRFDA call with per-variable fields
// Grid pointer must be passed (allocated via wrfda_alloc_and_init_domain)
// Observation location metadata is handled via iv structure
// Family is determined automatically from iv structure
// Output is stored internally and extracted via
// wrfda_extract_tangent_linear_output
int wrfda_xtoy_apply_grid(const void* grid_ptr, const void* ob_ptr,
                          const void* iv_ptr);

// Get count of tangent linear output values (count-only call)
int wrfda_get_tangent_linear_count(int* num_innovations);

// Extract output values from the tangent linear operator
int wrfda_extract_tangent_linear_output(double* out_y, int* num_innovations);

int wrfda_xtoy_adjoint_grid(const void* grid_ptr, const void* iv_ptr);

// Set delta_y input for adjoint operator
int wrfda_set_delta_y(const double* delta_y, int num_obs);

// Get adjoint gradients from persistent arrays
int wrfda_get_adjoint_gradients(double* u, double* v, double* t, double* q,
                                double* psfc);

// Enhanced function for constructing iv_type with detailed observation
// information
int wrfda_construct_iv_from_observations(
    const char* operator_family, int nx, int ny, int nz, const double* lats2d,
    const double* lons2d, const double* levels,  // size nx*ny, nz
    const int* num_obs, const double* obs_lats, const double* obs_lons,
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
int wrfda_get_innov_vector(const int* it, const void* ob_ptr,
                           const void* iv_ptr);

// Helper function to construct WRFDA domain structure from flat arrays
int wrfda_construct_domain_from_arrays(
    const int* nx, const int* ny, const int* nz, const double* u,
    const double* v, const double* t, const double* q, const double* psfc,
    const double* ph, const double* phb, const double* hf, const double* hgt,
    const double* p, const double* pb, const double* lats2d,
    const double* lons2d);

// Construct y_type from observation data
void* wrfda_construct_y_type(
    int* num_obs, int* num_levels, const double* u_values,
    const double* v_values, const double* t_values, const double* p_values,
    const double* q_values, const double* u_errors, const double* v_errors,
    const double* t_errors, const double* p_errors, const double* q_errors,
    const int* u_available, const int* v_available, const int* t_available,
    const int* p_available, const int* q_available, const double* lats,
    const double* lons, const char* obs_types, const char* family);

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
    const double* elevations, const char* obs_types, const char* family);

// Construct config_flags for WRFDA
void* wrfda_construct_config_flags();

// Count innovation values from iv_type structure
int wrfda_count_innovations(const char* family, int* num_innovations);

// Extract innovation values from iv_type structure
int wrfda_extract_innovations(const char* family, double* innovations,
                              int* num_innovations);

// Extract observation values from y_type structure
int wrfda_extract_observations(const char* family, double* observations,
                               int* num_observations);

// Initialize WRFDA variables for 3D-Var analysis
void initialize_wrfda_3dvar();

// Update analysis increments in WRFDA grid structure
void wrfda_update_analysis_increments(const double* u, const double* v,
                                      const double* t, const double* q,
                                      const double* psfc);

// Update background state (xb) from state data
void wrfda_update_background_state(const double* u, const double* v,
                                   const double* t, const double* q,
                                   const double* psfc, const double* ph,
                                   const double* phb, const double* hf,
                                   const double* hgt, const double* p,
                                   const double* pb, const double* lats2d,
                                   const double* lons2d);

// Get available observation families from iv structure
int wrfda_get_available_families(char* families_buffer, int* buffer_size);

// Cleanup da_control module vertical coordinates (should be called in
// destructor)
void wrfda_cleanup_vertical_coords();

#ifdef __cplusplus
}
#endif
