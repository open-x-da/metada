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

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Field information types (moved to wrfda_types.h for minimal dependencies)
#include "wrfda_types.h"

/**
 * @brief Tangent linear operator: H'(xb)·δx
 * @details Computes linearized observation operator for incremental 3D-Var.
 *          Applies WRFDA's da_transform_xtoy_* routines to grid%xa (increment).
 *          Output is extracted via wrfda_extract_observations for consistency.
 * @param grid_ptr Pointer to WRFDA grid structure (contains grid%xa)
 * @param ob_ptr Pointer to WRFDA y_type structure (output)
 * @param iv_ptr Pointer to WRFDA iv_type structure (observation metadata)
 * @return 0 on success, non-zero on error
 */
int wrfda_xtoy_apply_grid();

// Compute weighted residual and observation cost using WRFDA's proven workflow:
// 1. da_calculate_residual: re = (O-B) - H(xa)
// 2. da_calculate_grady: jo_grad_y = -R^{-1} · re
// 3. Jo = -0.5 * Σ(re · jo_grad_y) = 0.5 * re^T · R^{-1} · re
// Stores jo_grad_y in persistent_jo_grad_y for use by wrfda_xtoy_adjoint_grid
// Returns observation cost via jo_cost output parameter
int wrfda_compute_weighted_residual(const void* iv_ptr, const void* y_ptr,
                                    double* jo_cost);

// Compute ONLY observation cost without modifying persistent_jo_grad_y
// Used for cost-only evaluations during optimization
// Computes Jo = 0.5 * re^T * R^{-1} * re using temporary structures
int wrfda_compute_jo_only(const void* iv_ptr, const void* y_ptr,
                          double* jo_cost);

// Compute ONLY jo_grad_y without computing cost (for gradient-only operations)
// Used during CG iterations when only gradient is needed, not cost.
// Computes jo_grad_y = -R^{-1} · re where re = (O-B) - H'(δx)
// Stores result in persistent_jo_grad_y for use by adjoint operator
// This avoids expensive cost computation during minimization iterations
int wrfda_compute_grady_only(const void* iv_ptr, const void* y_ptr);

int wrfda_xtoy_adjoint_grid(const void* grid_ptr, const void* iv_ptr);

// Pure adjoint operator for TL/AD testing (without R^{-1} weighting)
// Takes raw dy values and applies H'^T directly
// Used for gradient checks: <H'dx, dy> = <dx, H'^T dy>
int wrfda_xtoy_adjoint_pure(const void* grid_ptr, const void* iv_ptr,
                            const double* dy_values, int num_values);

// Control-space adjoint entrypoint (if available in WRFDA build). Returns 0 on
// success; non-zero if not available or failed. Callers should fall back to
// state-adjoint + U^T when non-zero is returned.
int wrfda_vtoy_adjoint_control(const void* grid_ptr, const void* iv_ptr,
                               double* control_gradient_out, int control_size);

// Control-to-observation-space transform: y = H'(U·v)
// This is the forward transform used during CG iterations when computing
// gradients along search directions. Matches WRFDA's da_transform_vtoy.
// @param grid_ptr Pointer to WRFDA grid structure
// @param iv_ptr Pointer to innovation vector structure
// @param control_input Control variable vector v
// @param control_size Size of control vector
// @param y_ptr Pointer to y_type structure (output: y = H'(U·v))
// @return 0 on success, non-zero on error
int wrfda_transform_vtoy(const void* grid_ptr, const void* iv_ptr,
                         const double* control_input, int control_size,
                         void* y_ptr);

// (deprecated note) When wrfda_vtoy_adjoint_control is unavailable, use
// state-adjoint + U^T path instead.

// Get current WRFDA CV5 control vector size
int wrfda_control_backend_get_cv_size(void);

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
                           const void* iv_ptr, const void* grid_ptr);

// Helper function to construct WRFDA domain structure from flat arrays
int wrfda_construct_domain_from_arrays(
    const int* nx, const int* ny, const int* nz, const double* u,
    const double* v, const double* t, const double* q, const double* psfc,
    const double* ph, const double* phb, const double* hf, const double* hgt,
    const double* p, const double* pb, const double* lats2d,
    const double* lons2d);

// Count innovation values from iv_type structure for all observation types
int wrfda_count_innovations(const void* iv_ptr, int* num_innovations);

// Extract innovation values from iv_type structure for all families
int wrfda_extract_innovations(const void* iv_ptr, double* innovations,
                              int* num_innovations);

// Extract observation values from y_type structure for all families
int wrfda_extract_observations(const void* iv_ptr, const void* ob_ptr,
                               double* observations, int* num_observations);

// Initialize WRFDA variables for 3D-Var analysis
void initialize_wrfda_3dvar();

// Trace control helpers
bool wrfda_trace_is_enabled(void);
void wrfda_trace_initialize(void);
void wrfda_trace_finalize(void);

// Update all analysis increment fields in WRFDA grid%xa structure
// Includes 16 3D fields and 2 2D fields
void wrfda_update_analysis_increments(
    const double* u, const double* v, const double* w, const double* t,
    const double* p, const double* q, const double* qt, const double* rh,
    const double* rho, const double* geoh, const double* wh, const double* qcw,
    const double* qrn, const double* qci, const double* qsn, const double* qgr,
    const double* psfc, const double* mu, const void* grid_ptr);

/**
 * @brief Transfer WRF fields to background state structure (grid%xb)
 *
 * @details Wraps WRFDA's da_transfer_wrftoxb routine which:
 * - Transfers WRF native fields (u, v, t, q, etc.) to grid%xb structure
 * - Computes derived fields (pressure, height, etc.)
 * - Applies coordinate transformations for Arakawa-C grid
 * - Prepares grid%xb for use by observation operators
 *
 * This MUST be called before using any WRFDA observation operators to ensure
 * grid%xb is properly populated from the current state. This follows the
 * standard WRFDA workflow where da_transfer_wrftoxb is called before
 * da_get_innov_vector.
 *
 * @return 0 on success, non-zero on error
 *
 * @see da_transfer_wrftoxb.inc in WRFDA source
 * @see da_setup_firstguess_wrf.inc which calls da_transfer_wrftoxb
 */
int wrfda_transfer_wrftoxb(void);

// WARNING: DO NOT call in destructors!
// Cleanup da_control module vertical coordinates (module-level, shared by all
// instances) Should ONLY be called during final WRFDA shutdown, not
// per-instance cleanup
void wrfda_cleanup_vertical_coords();

// ============================================================================
// Observation Reading and Setup Functions
// ============================================================================

/**
 * @brief Read and allocate observations using WRFDA's standard pipeline
 * @details Reads observations from BUFR file and allocates WRFDA structures.
 *          Calls: da_setup_obs_structures_bufr → da_read_obs_bufr
 *          Does NOT compute innovations (happens later in ObsOperator)
 * @param[in] grid_ptr Pointer to WRFDA grid structure
 * @param[in] ob_filename Path to BUFR observation file
 * @param[in] ob_filename_len Length of filename string
 * @param[out] iv_ptr Output pointer to iv_type structure
 * @param[out] ob_ptr Output pointer to y_type structure
 * @return 0 on success, non-zero on error
 */
int wrfda_read_and_allocate_observations(void* grid_ptr,
                                         const char* ob_filename,
                                         int ob_filename_len, void** iv_ptr,
                                         void** ob_ptr);

/**
 * @brief Release observation structures allocated by
 * wrfda_read_and_allocate_observations.
 * @return 0 on success, non-zero on error.
 */
int wrfda_release_observations(void);

/**
 * @brief Get observation counts for all observation types
 * @param[in] iv_ptr Pointer to WRFDA iv_type structure
 * @param[out] obs_counts Array to receive counts (size: num_ob_indexes)
 * @return 0 on success, non-zero on error
 */
int wrfda_get_obs_type_counts(void* iv_ptr, int* obs_counts);

/**
 * @brief Get total observation count across all types
 * @param[in] iv_ptr Pointer to WRFDA iv_type structure
 * @param[out] total_count Total number of observations
 * @return 0 on success, non-zero on error
 */
int wrfda_get_total_obs_count(void* iv_ptr, int* total_count);

/**
 * @brief Get pointer to WRFDA iv_type structure (innovation vector)
 * @details Returns the C pointer to the module-level wrfda_iv structure
 *          that was allocated during wrfda_read_and_allocate_observations.
 *          This allows direct access to WRFDA's observation structures without
 *          going through backend objects.
 * @return C pointer to iv_type structure, or NULL if not allocated
 */
void* wrfda_get_iv_ptr(void);

/**
 * @brief Get pointer to WRFDA y_type structure (observation values)
 * @details Returns the C pointer to the module-level wrfda_ob structure
 *          that was allocated during wrfda_read_and_allocate_observations.
 *          This allows direct access to WRFDA's observation structures without
 *          going through backend objects.
 * @return C pointer to y_type structure, or NULL if not allocated
 */
void* wrfda_get_y_ptr(void);

#ifdef __cplusplus
}
#endif
