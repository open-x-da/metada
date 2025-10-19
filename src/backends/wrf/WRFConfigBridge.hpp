/**
 * @file WRFConfigBridge.hpp
 * @brief C++ interface to WRF Fortran configuration system
 * @ingroup backends
 * @author Metada Framework Team
 */

#pragma once

#include <cstdint>

extern "C" {

//============================================================================
// WRF Initialization Functions
//============================================================================

/**
 * @brief Initialize WRFU (WRF ESMF time utilities)
 *
 * @details Initializes the WRFU time manager with Gregorian calendar.
 * Must be called between wrf_init_modules_(1) and wrf_init_modules_(2).
 */
void wrf_wrfu_initialize_();

/**
 * @brief Initialize WRF modules
 *
 * @param[in] phase Initialization phase (1 or 2)
 * @details Phase 1 initializes core modules, Phase 2 initializes advanced
 * features
 */
void wrf_init_modules_(int phase);

/**
 * @brief Check if WRF has been initialized
 *
 * @return bool True if WRF modules have been initialized
 */
bool wrf_is_initialized_();

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
// WRF Domain Management Functions
//============================================================================

/**
 * @brief Allocate and initialize a WRF domain
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return int Error code (0 = success, non-zero = error)
 */
int wrf_alloc_domain_(int domain_id);

/**
 * @brief Deallocate a WRF domain
 *
 * @param[in] domain_id Domain ID (1-based)
 */
void wrf_dealloc_domain_(int domain_id);

/**
 * @brief Check if a domain exists
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return bool True if domain exists
 */
bool wrf_domain_exists_(int domain_id);

//============================================================================
// WRF Grid Geometry Setters
//============================================================================

/**
 * @brief Set grid geometry parameters in WRF domain structures
 *
 * @param[in] domain_id Domain ID (1-based)
 * @param[in] dx Grid spacing in X direction (meters)
 * @param[in] dy Grid spacing in Y direction (meters)
 * @param[in] map_proj Map projection type
 * @param[in] cen_lat Center latitude (degrees)
 * @param[in] cen_lon Center longitude (degrees)
 * @param[in] truelat1 First true latitude (degrees)
 * @param[in] truelat2 Second true latitude (degrees)
 * @param[in] stand_lon Standard longitude (degrees)
 */
void wrf_set_grid_geometry_(int domain_id, float dx, float dy, int map_proj,
                            float cen_lat, float cen_lon, float truelat1,
                            float truelat2, float stand_lon);

//============================================================================
// WRF Grid Geometry Getters
//============================================================================

/**
 * @brief Get dx from WRF domain structure
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Grid spacing in X direction (meters)
 */
float wrf_get_grid_dx_(int domain_id);

/**
 * @brief Get dy from WRF domain structure
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Grid spacing in Y direction (meters)
 */
float wrf_get_grid_dy_(int domain_id);

/**
 * @brief Get map_proj from WRF domain structure
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return int Map projection type
 */
int wrf_get_grid_map_proj_(int domain_id);

/**
 * @brief Get center latitude from WRF config_flags
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Center latitude (degrees)
 */
float wrf_get_grid_cen_lat_(int domain_id);

/**
 * @brief Get center longitude from WRF config_flags
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Center longitude (degrees)
 */
float wrf_get_grid_cen_lon_(int domain_id);

/**
 * @brief Get first true latitude from WRF config_flags
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float First true latitude (degrees)
 */
float wrf_get_grid_truelat1_(int domain_id);

/**
 * @brief Get second true latitude from WRF config_flags
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return float Second true latitude (degrees)
 */
float wrf_get_grid_truelat2_(int domain_id);

/**
 * @brief Get standard longitude from WRF config_flags
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

/**
 * @brief Get C pointer to WRF grid/domain structure
 *
 * @param[in] domain_id Domain ID (1-based)
 * @return void* Opaque pointer to domain structure, or nullptr if not found
 * @details This pointer can be passed to WRF observation operators and
 * other WRF Fortran routines that require grid as a parameter.
 *
 * @note The pointer is valid as long as the domain exists
 * @note Do not delete or free this pointer - it's managed by Fortran
 */
void* wrf_get_grid_ptr_(int domain_id);

}  // extern "C"

namespace metada::backends::wrf {

/**
 * @brief RAII wrapper for WRF configuration and domain initialization
 *
 * @details Manages complete WRF domain lifecycle including:
 * - WRF module initialization
 * - Namelist configuration reading
 * - Domain structure allocation
 * - Domain deallocation on destruction
 *
 * This ensures WRF's internal state is properly initialized and cleaned up
 * for WRFDA operations.
 *
 * @note This class does NOT provide map projection parameters. Those should
 * be obtained from WRFGeometry, which reads the authoritative values from
 * NetCDF file global attributes. The namelist.input values may be incorrect
 * or zeros, so always use WRFGeometry for grid parameters.
 *
 * @see WRFGeometry::getDx(), WRFGeometry::getMapProj(), etc.
 */
class WRFConfigManager {
 public:
  /**
   * @brief Initialize WRF configuration system and allocate domain
   *
   * @param[in] domain_id Domain ID (default: 1 for outermost domain)
   * @param[in] allocate_domain Whether to allocate WRF domain structure
   * (default: true)
   *
   * @details Performs the following steps:
   * 1. Initializes WRF modules (if not already initialized)
   * 2. Calls initial_config() to read namelist.input
   * 3. Calls model_to_grid_config_rec() for domain-specific config
   * 4. Optionally allocates WRF domain structure
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
   * @brief Destructor - deallocates WRF domain if allocated
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
   * @return bool True if WRF domain structure was allocated
   */
  bool isDomainAllocated() const { return domain_allocated_; }

  /**
   * @brief Check if domain exists in WRF
   *
   * @return bool True if WRF domain exists
   */
  bool domainExists() const;

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
   * @brief Get pointer to WRF grid/domain structure
   *
   * @return void* Opaque pointer to domain structure, or nullptr if not found
   * @details Returns a pointer to the WRF domain structure for this domain_id.
   * This can be passed to WRF observation operators.
   *
   * @note The pointer is valid as long as the domain exists
   * @note Do not delete or free this pointer - it's managed by Fortran
   */
  void* getGridPtr() const;

  /**
   * @note Map projection parameters (dx, dy, map_proj, etc.) are NOT provided
   * by this class. Those values should be obtained from WRFGeometry, which
   * reads them from the NetCDF file global attributes (the authoritative
   * source). This class only handles WRF internal configuration initialization.
   */

 private:
  int domain_id_;                  ///< Domain ID for this configuration
  bool initialized_ = false;       ///< Initialization status
  bool domain_allocated_ = false;  ///< Whether WRF domain was allocated

  /// Initialize WRF modules (called once globally)
  static void initializeWRFModules();

  /// Static flag to track global WRF initialization
  static bool wrf_modules_initialized_;
};

}  // namespace metada::backends::wrf
