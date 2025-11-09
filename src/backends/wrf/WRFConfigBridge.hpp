/**
 * @file WRFConfigBridge.hpp
 * @brief C++ interface to WRF Fortran configuration system
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <cstddef>

extern "C" {

//============================================================================
// WRFDA Initialization Functions
//============================================================================

/**
 * @brief Initialize WRFU (WRF ESMF time utilities)
 *
 * @details Initializes the WRFU time manager with Gregorian calendar.
 * Must be called between wrfda_init_modules_(1) and wrfda_init_modules_(2).
 */
void wrfda_wrfu_initialize_();

/**
 * @brief Initialize WRFDA modules
 *
 * @param[in] phase Initialization phase (1 or 2)
 * @details Phase 1 initializes core modules (configuration and constants),
 * Phase 2 initializes advanced features. Domain management is handled
 * separately by WRFDA itself.
 */
void wrfda_init_modules_(int phase);

/**
 * @brief Check if WRFDA has been initialized
 *
 * @return bool True if WRFDA modules have been initialized
 */
bool wrfda_is_initialized_();

/**
 * @brief Initialize WRF configuration by reading namelist.input
 *
 * @details This function calls WRF's initial_config() subroutine which:
 * - Opens and reads namelist.input file
 * - Populates the module-level model_config_rec structure
 * - Validates configuration parameters
 *
 * This must be called before any WRF operations that require configuration.
 *
 * @note The namelist.input file must exist in the current working directory
 * @note This populates the module-level model_config_rec in module_configure
 */
void wrf_initial_config_();

/**
 * @brief Copy namelist configuration from model_config_rec to da_control module
 *
 * @details This function performs the same configuration copy as WRFDA's
 * da_wrfvar_init1.inc by including config_assigns.inc. It copies all namelist
 * variables (like use_synopobs, use_metarobs, etc.) from model_config_rec
 * to the module-level variables in da_control.
 *
 * This MUST be called after wrf_initial_config_() and before any WRFDA
 * operations that access da_control module variables.
 *
 * @pre wrf_initial_config_() must have been called first
 */
void copy_config_to_da_control();

/**
 * @brief Validate WRFDA configuration for common conflicts and errors
 *
 * @details This function performs the same sanity checks as WRFDA's
 * da_solve.inc to catch configuration errors early. It validates:
 * - GPS observation conflicts (use_gpsrefobs vs use_gpsephobs)
 * - Radar observation conflicts (use_radar_rf vs use_radar_rhv)
 * - Control variable option validity
 * - Vertical interpolation parameters
 * - Cloud CV and alpha CV compatibility
 * - Dual-resolution hybrid constraints
 *
 * This catches user errors before they cause cryptic failures during DA.
 *
 * @param[out] error_code Error code (0=success, 1-9=specific validation
 * failure)
 *
 * @pre copy_config_to_da_control() must have been called first
 * @see da_solve.inc lines 168-252 in WRFDA source
 */
void validate_wrfda_config(int* error_code);

/**
 * @brief Query whether tracing is enabled via namelist.trace_use
 * @return bool True when trace_use=.true. in namelist configuration
 */
bool wrfda_trace_is_enabled(void);

/**
 * @brief Initialize WRFDA tracing subsystem if enabled
 */
void wrfda_trace_initialize(void);

/**
 * @brief Finalize WRFDA tracing subsystem and write reports
 */
void wrfda_trace_finalize(void);

/**
 * @brief Convert model config to domain-specific grid config
 *
 * @details This function calls WRF's model_to_grid_config_rec() subroutine
 * which extracts configuration for a specific domain from the model-wide
 * configuration and stores it in the module-level config_flags_ variable.
 *
 * @param[in] domain_id Domain ID (1-based, typically 1 for outermost domain)
 *
 * @pre wrf_initial_config_() must have been called first
 * @post Module-level config_flags_ contains domain-specific configuration
 *
 * @note This function is called during initialization for WRF internal
 * consistency. Map projection parameters should be obtained from WRFGeometry
 * which reads them from NetCDF file global attributes (the authoritative
 * source).
 */
void wrf_model_to_grid_config_(const int* domain_id);

//============================================================================
// WRFDA Domain Allocation and Management
//============================================================================

/**
 * @brief Allocate and initialize WRFDA domain following standard workflow
 *
 * @param[in] domain_id Domain ID (typically 1 for single-domain)
 * @return int Error code (0 = success, non-zero = error)
 *
 * @details Follows the exact sequence from da_wrfvar_init2.inc:
 * 1. alloc_and_configure_domain (allocates head_grid)
 * 2. model_to_grid_config_rec (extracts domain config)
 * 3. set_scalar_indices_from_config (sets tracer indices)
 * 4. init_wrfio (initializes WRF I/O system)
 * 5. setup_timekeeping (sets up domain clock)
 *
 * @note This must be called after wrf_initial_config_()
 */
int wrfda_alloc_and_init_domain_(int domain_id);

/**
 * @brief Get pointer to WRFDA's head_grid structure
 *
 * @return void* Pointer to head_grid, or nullptr if not allocated
 *
 * @note The pointer is valid as long as head_grid exists
 * @note Do not free this pointer - it's managed by WRFDA
 */
void* wrfda_get_head_grid_ptr_();

/**
 * @brief Check if WRFDA's head_grid has been allocated
 *
 * @return bool True if head_grid is allocated
 */
bool wrfda_head_grid_allocated_();

//============================================================================
// WRFDA Grid Geometry Getters
//============================================================================

/**
 * @brief Get dx from WRFDA's model_config_rec
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Grid spacing in X direction (meters)
 */
float wrf_get_grid_dx_(int domain_id);

/**
 * @brief Get dy from WRFDA's model_config_rec
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Grid spacing in Y direction (meters)
 */
float wrf_get_grid_dy_(int domain_id);

/**
 * @brief Get map_proj from WRFDA's model_config_rec
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return int Map projection type
 */
int wrf_get_grid_map_proj_(int domain_id);

/**
 * @brief Get center latitude from WRFDA's model_config_rec
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Center latitude (degrees)
 */
float wrf_get_grid_cen_lat_(int domain_id);

/**
 * @brief Get center longitude from WRFDA's model_config_rec
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Center longitude (degrees)
 */
float wrf_get_grid_cen_lon_(int domain_id);

/**
 * @brief Get first true latitude from WRFDA's model_config_rec
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float First true latitude (degrees)
 */
float wrf_get_grid_truelat1_(int domain_id);

/**
 * @brief Get second true latitude from WRFDA's model_config_rec
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Second true latitude (degrees)
 */
float wrf_get_grid_truelat2_(int domain_id);

/**
 * @brief Get standard longitude from WRFDA's model_config_rec
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Standard longitude (degrees)
 */
float wrf_get_grid_stand_lon_(int domain_id);

//============================================================================
// Config Flags Pointer Access
//============================================================================

/**
 * @brief Get C pointer to WRF config_flags structure
 *
 * @return void* Opaque pointer to grid_config_rec_type structure
 * @details This pointer can be passed to WRF observation operators and
 * other WRF Fortran routines that require config_flags as a parameter.
 *
 * @note The pointer points to module-level config_flags_ in
 * metada_wrf_config_bridge
 * @note The structure is maintained in Fortran - do not free this pointer
 */
void* wrf_get_config_flags_ptr_();

/**
 * @brief Get size of config_flags structure
 *
 * @return size_t Size of grid_config_rec_type structure in bytes
 * @details Useful for validation and debugging
 */
size_t wrf_get_config_flags_size_();

//============================================================================
// State Operations - Zero Increment
//============================================================================

/**
 * @brief Zero the analysis increment (xa) using WRFDA's da_zero_x
 *
 * @details Wraps WRFDA's proven da_zero_x subroutine which zeros all fields
 * in the analysis increment structure (grid%xa). This includes:
 * - 3D fields: u, v, w, t, q, p, geoh, rh, wh, rho, ref
 * - Hydrometeors: qcw, qrn, qt, qci, qsn, qgr
 * - 2D fields: tgrn, psfc, mu, u10, v10, t2, q2
 * - Derived fields: ztd, tpw, speed, brightness temperatures
 *
 * @param[in] grid_ptr Pointer to WRF domain structure
 * @return int Error code: 0 on success, -1 if grid not initialized
 *
 * @note The grid must be initialized via wrfda_alloc_and_init_domain
 *
 * @see da_zero_x in WRF/var/da/da_define_structures/da_zero_x.inc
 */
int wrfda_zero_xa(void* grid_ptr);

}  // extern "C"

namespace metada::backends::wrf {

/**
 * @brief RAII wrapper for WRFDA initialization and domain allocation
 *
 * @details Manages complete WRFDA initialization lifecycle following the
 * standard workflow from da_wrfvar_init1.inc and da_wrfvar_init2.inc:
 * - WRFDA module initialization
 * - Namelist configuration reading
 * - Domain allocation and configuration
 * - WRF I/O initialization
 * - Timekeeping setup
 *
 * This ensures WRFDA's internal state is properly initialized for
 * observation operators and data assimilation operations.
 *
 * @note Grid geometry parameters come from namelist.input, which MUST match
 * the NetCDF input file. Use WRFGeometry to verify/validate these values.
 *
 * @see WRFGeometry for grid parameter validation
 */
class WRFConfigManager {
 public:
  /**
   * @brief Initialize WRFDA configuration system and allocate domain
   *
   * @param[in] domain_id Domain ID (default: 1 for outermost domain)
   * @param[in] allocate_domain Whether to allocate domain (default: true)
   *
   * @details Performs the following steps (matching WRFDA workflow):
   * 1. Initializes WRFDA modules (if not already initialized)
   * 2. Calls initial_config() to read namelist.input
   * 3. If allocate_domain:
   *    a. Calls alloc_and_configure_domain()
   *    b. Calls model_to_grid_config_rec()
   *    c. Calls set_scalar_indices_from_config()
   *    d. Calls init_wrfio()
   *    e. Calls setup_timekeeping()
   *
   * @throws std::runtime_error If initialization fails
   */
  explicit WRFConfigManager(int domain_id = 1, bool allocate_domain = true);

  // Prevent copying
  WRFConfigManager(const WRFConfigManager&) = delete;
  WRFConfigManager& operator=(const WRFConfigManager&) = delete;

  // Allow moving
  WRFConfigManager(WRFConfigManager&&) noexcept = default;
  WRFConfigManager& operator=(WRFConfigManager&&) noexcept = default;

  /**
   * @brief Destructor - cleanup (no domain deallocation needed)
   *
   * @note WRFDA's head_grid persists and is not deallocated
   */
  ~WRFConfigManager();

  /**
   * @brief Get domain ID
   *
   * @return int Domain ID (1-based)
   */
  int getDomainId() const { return domain_id_; }

  /**
   * @brief Check if domain was allocated
   *
   * @return bool True if WRFDA domain structure was allocated
   */
  bool isDomainAllocated() const { return domain_allocated_; }

  /**
   * @brief Get pointer to WRFDA's head_grid structure
   *
   * @return void* Pointer to head_grid, or nullptr if not allocated
   * @details This pointer can be passed to WRFDA observation operators
   * and other WRFDA Fortran routines that require grid as a parameter.
   *
   * @note The pointer is valid as long as the domain exists
   * @note Do not delete or free this pointer - it's managed by WRFDA
   */
  void* getGridPtr() const;

  /**
   * @brief Get pointer to WRF config_flags structure
   *
   * @return void* Opaque pointer to grid_config_rec_type structure
   * @details This pointer can be passed to WRF observation operators and
   * other WRF Fortran routines that require config_flags as a parameter.
   *
   * @note The pointer is valid as long as the WRFConfigManager instance exists
   * @note Do not delete or free this pointer - it's managed by Fortran
   *
   * Example usage:
   * @code
   * WRFConfigManager config(1);
   * void* config_flags_ptr = config.getConfigFlagsPtr();
   * // Pass to WRF observation operator:
   * call_wrf_obsoperator(grid_ptr, config_flags_ptr, ...);
   * @endcode
   */
  void* getConfigFlagsPtr() const;

  /**
   * @brief Get size of config_flags structure
   *
   * @return size_t Size in bytes
   * @details Useful for validation and debugging
   */
  size_t getConfigFlagsSize() const;

  /**
   * @note Map projection parameters (dx, dy, map_proj, etc.) are NOT provided
   * by this class. Those values should be obtained from WRFGeometry, which
   * reads them from the NetCDF file global attributes (the authoritative
   * source). This class only handles WRFDA internal initialization.
   */

 private:
  int domain_id_;                  ///< Domain ID for this configuration
  bool initialized_ = false;       ///< Initialization status
  bool domain_allocated_ = false;  ///< Whether WRFDA domain was allocated
  bool trace_session_started_ =
      false;  ///< Whether tracing was initialized for this session

  /// Initialize WRFDA modules (called once globally)
  static void initializeWRFDAModules();

  /// Static flag to track global WRFDA initialization
  static bool wrfda_modules_initialized_;
};

}  // namespace metada::backends::wrf
