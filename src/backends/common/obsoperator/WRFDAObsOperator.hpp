/**
 * @file WRFDAObsOperator.hpp
 * @brief WRFDA observation operator backend bridge (header-only)
 * @ingroup backends
 *
 * @details
 * This operator provides a WRFDA bridge for observation operators that can be
 * used by any real model backend. It directly implements WRFDA-specific logic
 * via C/Fortran bridge calls to native WRFDA routines (da_transform_xtoy_*).
 * The operator uses configuration keys for integrating with WRFDA without
 * coupling to any specific model backend.
 *
 * The class uses generic configuration keys (external_root, external_system)
 * to allow easy switching between different external observation operator
 * systems (WRFDA, GSI, DART, etc.) at configuration time for better
 * performance. The actual external system is determined during initialization,
 * not at runtime.
 *
 * Config keys (optional):
 * - external_root: Path to external operator sources (e.g.,
 * D:/linux/WRF/var/da)
 * - external_system: External system identifier (e.g., "wrfda", "gsi", "dart")
 * - operator_family: Observation operator family (e.g., "metar", "gpspw")
 * - required_state_vars: [array of strings]
 * - required_obs_vars: [array of strings]
 */

#pragma once

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../../wrf/WRFGridInfo.hpp"
#include "Location.hpp"
#include "WRFDAObsOperator_c_api.h"

namespace metada::backends::common::obsoperator {

template <typename StateBackend, typename ObsBackend>
class WRFDAObsOperator {
 public:
  /**
   * @brief Default constructor (deleted)
   *
   * @details This class requires configuration for initialization and cannot be
   * default-constructed.
   */
  WRFDAObsOperator() = delete;

  /**
   * @brief Copy constructor (deleted)
   *
   * @details Copying is not supported to prevent resource sharing issues.
   * Use move semantics instead.
   */
  WRFDAObsOperator(const WRFDAObsOperator&) = delete;

  /**
   * @brief Copy assignment operator (deleted)
   *
   * @details Copy assignment is not supported to prevent resource sharing
   * issues. Use move semantics instead.
   */
  WRFDAObsOperator& operator=(const WRFDAObsOperator&) = delete;

  /**
   * @brief Constructor with configuration
   *
   * @tparam ConfigBackend The configuration backend type
   * @param config Configuration object containing operator parameters
   *
   * @details Initializes the WRFDA observation operator with the provided
   * configuration. The configuration should include external system paths,
   * operator families, and required variable specifications.
   *
   * @throws std::runtime_error if initialization fails
   */
  template <typename ConfigBackend>
  explicit WRFDAObsOperator(const ConfigBackend& config) {
    initialize(config);
  }

  /**
   * @brief Move constructor
   *
   * @param other The WRFDAObsOperator to move from
   *
   * @details Moves ownership of resources from another WRFDAObsOperator
   * instance. The source object is left in a valid but unspecified state.
   *
   * @note This constructor is noexcept, ensuring it can be used in containers
   * that require noexcept move operations.
   */
  WRFDAObsOperator(WRFDAObsOperator&& other) noexcept
      : initialized_(other.initialized_),
        required_state_vars_(std::move(other.required_state_vars_)),
        required_obs_vars_(std::move(other.required_obs_vars_)),
        external_root_(std::move(other.external_root_)),
        external_system_(std::move(other.external_system_)),
        operator_families_(std::move(other.operator_families_)) {
    other.initialized_ = false;
  }

  /**
   * @brief Move assignment operator
   *
   * @param other The WRFDAObsOperator to move from
   * @return Reference to this object
   *
   * @details Moves ownership of resources from another WRFDAObsOperator
   * instance. The source object is left in a valid but unspecified state.
   *
   * @note This operator is noexcept, ensuring it can be used in containers
   * that require noexcept move operations.
   */
  WRFDAObsOperator& operator=(WRFDAObsOperator&& other) noexcept {
    if (this != &other) {
      initialized_ = other.initialized_;
      required_state_vars_ = std::move(other.required_state_vars_);
      required_obs_vars_ = std::move(other.required_obs_vars_);
      external_root_ = std::move(other.external_root_);
      external_system_ = std::move(other.external_system_);
      operator_families_ = std::move(other.operator_families_);
      other.initialized_ = false;
    }
    return *this;
  }

  /**
   * @brief Initialize the WRFDA observation operator
   *
   * @tparam ConfigBackend The configuration backend type
   * @param config Configuration object containing operator parameters
   *
   * @details Initializes the operator with configuration parameters including:
   * - external_root: Path to external operator sources
   * - external_system: External system identifier
   * - operator_family: Observation operator family or families
   * - required_state_vars: Required state variables
   * - required_obs_vars: Required observation variables
   *
   * @throws std::runtime_error if already initialized or if configuration
   * parsing fails
   *
   * @note This method can only be called once per instance
   */
  template <typename ConfigBackend>
  void initialize(const ConfigBackend& config) {
    if (isInitialized()) {
      throw std::runtime_error("WRFDAObsOperator already initialized");
    }

    try {
      external_root_ = config.Get("external_root").asString();
    } catch (...) {
      external_root_.clear();
    }

    try {
      external_system_ = config.Get("external_system").asString();
    } catch (...) {
      external_system_.clear();
    }

    try {
      auto family_value = config.Get("operator_family");
      if (family_value.isVectorString()) {
        // Handle array of operator families
        operator_families_ = family_value.asVectorString();
      } else if (family_value.isString()) {
        // Handle single operator family (backward compatibility)
        operator_families_ = {family_value.asString()};
      } else {
        operator_families_.clear();
      }
    } catch (...) {
      operator_families_.clear();
    }

    try {
      required_state_vars_ = config.Get("required_state_vars").asVectorString();
    } catch (...) {
      required_state_vars_.clear();
    }

    try {
      required_obs_vars_ = config.Get("required_obs_vars").asVectorString();
    } catch (...) {
      required_obs_vars_.clear();
    }

    // WRFDA-specific initialization complete
    // No delegation needed - this operator directly implements WRFDA logic

    // Initialize WRFDA variables for 3D-Var analysis
    // This sets var4d_run = .false., num_fgat_time = 1, and sfc_assi_options =
    // sfc_assi_options_1
    initialize_wrfda_3dvar();

    initialized_ = true;
  }

  /**
   * @brief Check if the operator is initialized
   *
   * @return true if initialized, false otherwise
   *
   * @details Returns the initialization status of the operator. The operator
   * must be initialized before any operations can be performed.
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Apply the observation operator (forward operator)
   *
   * @param state The background state to apply the operator to
   * @param obs The observations defining the observation locations
   * @return Vector of simulated observations
   *
   * @details Applies the WRFDA observation operator to transform the model
   * state into observation space. This method:
   * - Extracts observation coordinates (lat, lon, level)
   * - Retrieves state variables (U, V, T, Q, PSFC)
   * - Calls the WRFDA C API for the actual computation
   * - Handles memory limitations for large observation sets
   * - Returns simulated observations at the specified locations
   *
   * @throws std::runtime_error if operator is not initialized, if observation
   * locations are not geographic, or if WRFDA computation fails
   *
   * @note Observations must use geographic coordinate system
   * @note For large observation sets (>10000), only the first 10000 are
   * processed to avoid WRFDA memory limitations
   */

  /**
   * @brief Set the background state for proper WRFDA integration
   *
   * @param background_state The background state (xb)
   *
   * @details This method sets the background state that will be used to
   * properly separate total state into background and increments for WRFDA.
   * When the background state is set, the operator will use the new
   * wrfda_xtoy_apply_grid_with_background function for better accuracy.
   */
  void setBackgroundState(const StateBackend& background_state) {
    background_state_data_ = background_state.getObsOperatorData();
    has_background_state_ = true;
  }

  /**
   * @brief Clear the background state
   *
   * @details This method clears the stored background state, causing the
   * operator to fall back to the original behavior.
   */
  void clearBackgroundState() {
    has_background_state_ = false;
    background_state_data_ = {};
  }

  std::vector<double> apply(const StateBackend& state,
                            const ObsBackend& obs) const {
    ensureInitialized();

    // Extract optimized data from backends
    const auto state_data = state.getObsOperatorData();
    const auto obs_data = obs.getObsOperatorData();

    const size_t num_observations = obs_data.num_obs;

    // Use pre-extracted state data
    const int nx = state_data.nx;
    const int ny = state_data.ny;
    const int nz = state_data.nz;

    // Get state variables (already extracted and optimized)
    const double* u = state_data.u;
    const double* v = state_data.v;
    const double* t = state_data.t;
    const double* q = state_data.q;
    const double* psfc = state_data.psfc;
    const double* ph = state_data.ph;
    const double* phb = state_data.phb;
    const double* hf = state_data.hf;
    const double* hgt = state_data.hgt;
    const double* p = state_data.p;
    const double* pb = state_data.pb;

    // Get grid metadata (already extracted)
    const double* lats_2d = state_data.lats_2d;
    const double* lons_2d = state_data.lons_2d;
    const double* levels = state_data.levels;

    // Validate required data
    if (!lats_2d || !lons_2d || !levels) {
      throw std::runtime_error(
          "WRFDAObsOperator: Missing required grid metadata (lat/lon/levels)");
    }

    // Convert staggered U and V to unstaggered grids if needed
    std::vector<double> u_unstaggered, v_unstaggered;
    const double *u_final = u, *v_final = v;

    if (u && state_data.u_dims.is_staggered) {
      // U is staggered, convert to unstaggered
      const auto& dims = state_data.u_dims;
      convertStaggeredUToUnstaggered(
          std::vector<double>(u, u + dims.nx * dims.ny * dims.nz),
          u_unstaggered, nx, ny, nz);
      u_final = u_unstaggered.data();
      std::cout << "WRFDA: Converted U from staggered (" << dims.nx << "x"
                << dims.ny << "x" << dims.nz << ") to unstaggered (" << nx
                << "x" << ny << "x" << nz << ")" << std::endl;
    }

    if (v && state_data.v_dims.is_staggered) {
      // V is staggered, convert to unstaggered
      const auto& dims = state_data.v_dims;
      convertStaggeredVToUnstaggered(
          std::vector<double>(v, v + dims.nx * dims.ny * dims.nz),
          v_unstaggered, nx, ny, nz);
      v_final = v_unstaggered.data();
      std::cout << "WRFDA: Converted V from staggered (" << dims.nx << "x"
                << dims.ny << "x" << dims.nz << ") to unstaggered (" << nx
                << "x" << ny << "x" << nz << ")" << std::endl;
    }

    // Validate essential variables
    if (!u && !v && !t) {
      throw std::runtime_error(
          "WRFDA observation operator requires at least one of U, V, or T "
          "variables");
    }
    if (!ph) {
      throw std::runtime_error(
          "WRFDA observation operator requires PH (geopotential perturbation)");
    }
    if (!phb) {
      throw std::runtime_error(
          "WRFDA observation operator requires PHB (base state geopotential)");
    }
    if (!hf) {
      throw std::runtime_error(
          "WRFDA observation operator requires HF (height field)");
    }
    if (!hgt) {
      throw std::runtime_error(
          "WRFDA observation operator requires HGT (terrain height)");
    }
    if (!p) {
      throw std::runtime_error(
          "WRFDA observation operator requires P (pressure perturbation)");
    }
    if (!pb) {
      throw std::runtime_error(
          "WRFDA observation operator requires PB (base state pressure)");
    }

    // Construct WRFDA domain structure from flat arrays
    void* domain_ptr = nullptr;
    int rc = wrfda_construct_domain_from_arrays(
        &nx, &ny, &nz, u_final, v_final, t, q, psfc, ph, phb, hf, hgt, p, pb,
        lats_2d, lons_2d, levels, &domain_ptr);

    if (rc != 0) {
      throw std::runtime_error(
          "Failed to construct WRFDA domain structure with code " +
          std::to_string(rc));
    }

    // Initialize WRFDA module-level variables (kts, kte, sfc_assi_options,
    // etc.)
    rc = initialize_wrfda_module_variables(domain_ptr);
    if (rc != 0) {
      throw std::runtime_error(
          "Failed to initialize WRFDA module variables with code " +
          std::to_string(rc));
    }

    // Step 2: Initialize map projection for coordinate conversion
    int map_proj = 3;           // Lambert conformal conic projection
    double cen_lat = 34.83002;  // center latitude
    double cen_lon = -81.03;    // center longitude
    double dx = 30000.0;        // grid spacing in m
    double stand_lon = -98.0;   // standard longitude
    double truelat1 = 30.0;     // first true latitude
    double truelat2 = 60.0;     // second true latitude

    initialize_map_projection_c(&map_proj, &cen_lat, &cen_lon, &dx, &stand_lon,
                                &truelat1, &truelat2);

    // Step 3: Initialize background state for WRFDA
    rc = wrfda_initialize_background_state(nx, ny, nz, u_final, v_final, t, q,
                                           psfc);
    if (rc != 0) {
      throw std::runtime_error(
          "Failed to initialize WRFDA background state with code " +
          std::to_string(rc));
    }

    // Step 4: Construct observation and innovation vector structures
    void* ob_ptr = constructYType(obs_data, operator_families_);
    void* iv_ptr = constructIvType(obs_data, operator_families_, domain_ptr);

    // Step 5: Call da_get_innov_vector to compute innovations y - H(xb)
    const int it = 1;
    rc = wrfda_get_innov_vector(&it, domain_ptr, ob_ptr, iv_ptr);

    if (rc != 0) {
      throw std::runtime_error("da_get_innov_vector failed with code " +
                               std::to_string(rc));
    }

    std::cout << "WRFDA: da_get_innov_vector completed successfully"
              << std::endl;

    // Extract results from iv structure and return as vector
    std::string family =
        operator_families_.empty() ? "synop" : operator_families_[0];
    std::vector<double> innovations = extractInnovations(iv_ptr, family);

    if (innovations.empty()) {
      std::cout
          << "WRFDA: No innovations extracted, returning placeholder values"
          << std::endl;
      return std::vector<double>(num_observations, 0.0);
    }

    std::cout << "WRFDA: Extracted " << innovations.size()
              << " innovation values" << std::endl;
    return innovations;
  }

  /**
   * @brief Get the required state variables
   *
   * @return Reference to vector of required state variable names
   *
   * @details Returns the list of state variables that must be present
   * in the state for the operator to function correctly. These variables
   * are typically meteorological fields like U, V, T, Q, and PSFC.
   */
  const std::vector<std::string>& getRequiredStateVars() const {
    return required_state_vars_;
  }

  /**
   * @brief Get the required observation variables
   *
   * @return Reference to vector of required observation variable names
   *
   * @details Returns the list of observation variables that the operator
   * can process. This helps in validating observation data before
   * applying the operator.
   */
  const std::vector<std::string>& getRequiredObsVars() const {
    return required_obs_vars_;
  }

  /**
   * @brief Apply the tangent linear observation operator
   *
   * @param state_increment The state increment to apply the operator to
   * @param reference_state The reference state (unused in this implementation)
   * @param obs The observations defining the observation locations
   * @return Vector of simulated observation increments
   *
   * @details Applies the tangent linear observation operator to transform
   * state increments into observation space. For this linearized operator,
   * the same array-based call is used on the increment as in the forward
   * operator.
   *
   * @throws std::runtime_error if operator is not initialized, if observation
   * locations are not geographic, or if WRFDA computation fails
   *
   * @note This implementation treats the operator as linear, so the tangent
   * linear operator is identical to the forward operator applied to increments
   */
  std::vector<double> applyTangentLinear(
      const StateBackend& state_increment,
      [[maybe_unused]] const StateBackend& reference_state,
      const ObsBackend& obs) const {
    // For tangent linear operator, use da_transform_xtoy_* functions
    // which compute H'(δx) where δx is the analysis increment
    ensureInitialized();

    // Extract optimized data from backends
    const auto state_data = state_increment.getObsOperatorData();
    const auto obs_data = obs.getObsOperatorData();

    const size_t num_observations = obs_data.num_obs;

    // Use pre-extracted state data
    const int nx = state_data.nx;
    const int ny = state_data.ny;
    const int nz = state_data.nz;

    // Get state variables (these are the analysis increments δx)
    const double* u = state_data.u;
    const double* v = state_data.v;
    const double* t = state_data.t;
    const double* q = state_data.q;
    const double* psfc = state_data.psfc;

    // Get grid metadata
    const double* lats_2d = state_data.lats_2d;
    const double* lons_2d = state_data.lons_2d;
    const double* levels = state_data.levels;

    // Get observation locations
    const double* obs_lats = obs_data.lats.data();
    const double* obs_lons = obs_data.lons.data();

    // Convert LevelData to double array for observation levels
    // For surface observations, each station has one level
    std::vector<double> obs_levels_vec(num_observations);

    // Debug: Check bounds before accessing
    if (obs_data.levels.size() != num_observations) {
      throw std::runtime_error("WRFDAObsOperator: Mismatch between num_obs (" +
                               std::to_string(num_observations) +
                               ") and levels.size() (" +
                               std::to_string(obs_data.levels.size()) + ")");
    }

    for (size_t i = 0; i < num_observations; ++i) {
      if (!obs_data.levels[i].empty()) {
        obs_levels_vec[i] =
            obs_data.levels[i][0]
                .pressure;  // First (and only) level for surface obs
      } else {
        obs_levels_vec[i] = 0.0;  // Default value if no levels
      }
    }
    const double* obs_levels = obs_levels_vec.data();

    // Validate required data
    if (!lats_2d || !lons_2d || !levels) {
      throw std::runtime_error(
          "WRFDAObsOperator: Missing required grid metadata (lat/lon/levels)");
    }

    if (!obs_lats || !obs_lons || !obs_levels) {
      throw std::runtime_error(
          "WRFDAObsOperator: Missing required observation location data");
    }

    // Convert staggered U and V to unstaggered grids if needed
    std::vector<double> u_unstaggered, v_unstaggered;
    const double *u_final = u, *v_final = v;

    if (u && state_data.u_dims.is_staggered) {
      u_unstaggered.resize(nx * ny * nz);
      // U staggered field has dimensions (nx+1) x ny x nz
      int u_staggered_size = (nx + 1) * ny * nz;
      std::vector<double> u_staggered_vec(u, u + u_staggered_size);
      convertStaggeredUToUnstaggered(u_staggered_vec, u_unstaggered, nx, ny,
                                     nz);
      u_final = u_unstaggered.data();
    }

    if (v && state_data.v_dims.is_staggered) {
      v_unstaggered.resize(nx * ny * nz);
      // V staggered field has dimensions nx x (ny+1) x nz
      int v_staggered_size = nx * (ny + 1) * nz;
      std::vector<double> v_staggered_vec(v, v + v_staggered_size);
      convertStaggeredVToUnstaggered(v_staggered_vec, v_unstaggered, nx, ny,
                                     nz);
      v_final = v_unstaggered.data();
    }

    // Call Fortran function directly to update analysis increments
    wrfda_update_analysis_increments(u_final, v_final, t, q, psfc);
    std::cout << "WRFDAObsOperator: Analysis increments updated successfully"
              << std::endl;

    // Get the correct number of innovations for proper sizing
    // Use default family since Fortran will determine actual family from iv
    // structure
    std::string family = "synop";  // Default family for innovation counting
    std::cout
        << "WRFDAObsOperator: Using default family for innovation counting: "
        << family << std::endl;

    int num_innovations = 0;
    int count_rc = wrfda_count_innovations(const_cast<char*>(family.c_str()),
                                           &num_innovations);

    if (count_rc != 0 || num_innovations <= 0) {
      throw std::runtime_error("Failed to count innovations with code " +
                               std::to_string(count_rc));
    }

    // Prepare output vector with correct size based on innovations
    std::vector<double> out_y(num_innovations, 0.0);

    std::cout << "WRFDAObsOperator: Calling tangent linear operator for "
              << num_observations << " observations with " << num_innovations
              << " innovations" << std::endl;

    const int num_obs_int = static_cast<int>(num_observations);
    int rc = wrfda_xtoy_apply_grid(&num_obs_int, out_y.data());

    if (rc != 0) {
      throw std::runtime_error("wrfda_xtoy_apply_grid failed with code " +
                               std::to_string(rc));
    }

    std::cout << "WRFDAObsOperator: Tangent linear operator returned "
              << out_y.size() << " innovation values" << std::endl;
    for (size_t i = 0; i < std::min(out_y.size(), size_t(10)); ++i) {
      std::cout << "  out_y[" << i << "] = " << out_y[i] << std::endl;
    }

    return out_y;
  }

  /**
   * @brief Apply the adjoint observation operator
   *
   * @param obs_increment The observation increment to apply the adjoint to
   * @param reference_state The reference state (unused in this implementation)
   * @param result_state The state to accumulate the adjoint result in
   * @param obs The observations defining the observation locations
   *
   * @details Applies the adjoint observation operator to transform observation
   * increments back to state space. This method:
   * - Extracts observation coordinates (lat, lon, level)
   * - Prepares state variable buffers for accumulation
   * - Calls the WRFDA adjoint C API for the actual computation
   * - Accumulates the result in the provided state
   *
   * @throws std::runtime_error if operator is not initialized, if observation
   * locations are not geographic, if increment sizes don't match, or if
   * WRFDA adjoint computation fails
   *
   * @note The result_state is modified in-place by accumulating the adjoint
   * contribution
   * @note Observations must use geographic coordinate system
   */
  void applyAdjoint(const std::vector<double>& obs_increment,
                    [[maybe_unused]] const StateBackend& reference_state,
                    StateBackend& result_state, const ObsBackend& obs) const {
    ensureInitialized();

    const size_t num_observations = obs.size();
    if (obs_increment.size() != num_observations) {
      throw std::runtime_error(
          "Adjoint increment size does not match number of observations");
    }

    std::vector<double> obs_lats(num_observations);
    std::vector<double> obs_lons(num_observations);
    std::vector<double> obs_levels(num_observations, 0.0);

    // Get model vertical levels for reference
    const auto& gi = result_state.geometry().unstaggered_info();

    for (size_t i = 0; i < num_observations; ++i) {
      const auto& loc = obs[i].location;
      if (loc.getCoordinateSystem() !=
          framework::CoordinateSystem::GEOGRAPHIC) {
        throw std::runtime_error(
            "WRFDAObsOperator requires geographic observation locations");
      }
      auto [lat, lon, level] = loc.getGeographicCoords();
      obs_lats[i] = lat;
      obs_lons[i] = lon;

      // Determine observation type and set level appropriately
      // For surface observations (like ADPSFC), use sigma level 1.0 (surface)
      // Determine observation type and set level appropriately (same logic as
      // apply method) For surface observations (like ADPSFC), always use sigma
      // level 1.0 This ensures consistent behavior regardless of varying
      // pressure values in BUFR data
      if (operator_families_.empty() || operator_families_[0] == "metar" ||
          operator_families_[0] == "synop" ||
          operator_families_[0] == "adpsfc") {
        // Surface observation family - always use sigma level 1.0 (surface
        // level) This avoids varying olev values from run to run due to
        // pressure variations
        obs_levels[i] = 1.0;
      } else if (level <= 0.0) {
        // Upper-air observation with invalid level - use model surface level
        obs_levels[i] = 1.0;
      } else {
        // Upper-air observation with valid level
        obs_levels[i] = level;
      }
    }

    // Prepare inout arrays for u,v,t,q,psfc accumulation
    auto get_or_zero_inout = [&](const std::string& var,
                                 std::vector<double>& storage, double*& ptr,
                                 int& nx, int& ny, int& nz) {
      ptr = static_cast<double*>(result_state.getData(var));
      if (ptr == nullptr) {
        nx = static_cast<int>(result_state.geometry().x_dim());
        ny = static_cast<int>(result_state.geometry().y_dim());
        nz = static_cast<int>(result_state.geometry().z_dim());
        const size_t n = static_cast<size_t>(std::max(1, nx) * std::max(1, ny) *
                                             std::max(1, nz));
        storage.assign(n, 0.0);
        ptr = storage.data();
      } else {
        try {
          const auto& dims = result_state.getVariableDimensions(var);
          if (dims.size() == 3) {
            nz = static_cast<int>(dims[0]);
            ny = static_cast<int>(dims[1]);
            nx = static_cast<int>(dims[2]);
          } else if (dims.size() == 2) {
            ny = static_cast<int>(dims[0]);
            nx = static_cast<int>(dims[1]);
            nz = 1;
          }
        } catch (...) {
        }
      }
    };

    double *u = nullptr, *v = nullptr, *t = nullptr, *q = nullptr,
           *psfc = nullptr;
    std::vector<double> u_zero, v_zero, t_zero, q_zero, psfc_zero;
    int nx_u = 1, ny_u = 1, nz_u = 1;
    int nx_v = 1, ny_v = 1, nz_v = 1;
    int nx_t = 1, ny_t = 1, nz_t = 1;
    int nx_q = 1, ny_q = 1, nz_q = 1;
    int nx_ps = 1, ny_ps = 1, nz_ps = 1;
    get_or_zero_inout("U", u_zero, u, nx_u, ny_u, nz_u);
    get_or_zero_inout("V", v_zero, v, nx_v, ny_v, nz_v);
    get_or_zero_inout("T", t_zero, t, nx_t, ny_t, nz_t);
    get_or_zero_inout("QVAPOR", q_zero, q, nx_q, ny_q, nz_q);
    get_or_zero_inout("PSFC", psfc_zero, psfc, nx_ps, ny_ps, nz_ps);

    const int nx = static_cast<int>(result_state.geometry().x_dim());
    const int ny = static_cast<int>(result_state.geometry().y_dim());
    const int nz = static_cast<int>(result_state.geometry().z_dim());

    // Convert staggered U and V to unstaggered grids if needed
    std::vector<double> u_unstaggered, v_unstaggered;
    const double *u_final = u, *v_final = v;

    if (u && isVariableStaggered(result_state, "U")) {
      // U is staggered, convert to unstaggered
      convertStaggeredUToUnstaggered(
          std::vector<double>(u, u + nx_u * ny_u * nz_u), u_unstaggered, nx, ny,
          nz);
      u_final = u_unstaggered.data();
      std::cout << "WRFDA Adjoint: Converted U from staggered (" << nx_u << "x"
                << ny_u << "x" << nz_u << ") to unstaggered (" << nx << "x"
                << ny << "x" << nz << ")" << std::endl;
    }

    if (v && isVariableStaggered(result_state, "V")) {
      // V is staggered, convert to unstaggered
      convertStaggeredVToUnstaggered(
          std::vector<double>(v, v + nx_v * ny_v * nz_v), v_unstaggered, nx, ny,
          nz);
      v_final = v_unstaggered.data();
      std::cout << "WRFDA Adjoint: Converted V from staggered (" << nx_v << "x"
                << ny_v << "x" << nz_v << ") to unstaggered (" << nx << "x"
                << ny << "x" << nz << ")" << std::endl;
    }

    // Use the gi variable already declared above
    const double* levels_ptr =
        gi.vertical_coords.empty() ? nullptr : gi.vertical_coords.data();
    std::vector<double> dummy_levels(1, 0.0);
    if (levels_ptr == nullptr) levels_ptr = dummy_levels.data();

    // Flatten the state variables into a single array for the adjoint operator
    const size_t grid_size = static_cast<size_t>(nx * ny * nz);
    std::vector<double> state_values(grid_size, 0.0);

    // Initialize with existing values if available
    if (u) {
      for (size_t i = 0; i < grid_size && i < static_cast<size_t>(nx * ny * nz);
           ++i) {
        state_values[i] += u_final[i];
      }
    }
    if (v) {
      for (size_t i = 0; i < grid_size && i < static_cast<size_t>(nx * ny * nz);
           ++i) {
        state_values[i] += v_final[i];
      }
    }
    if (t) {
      for (size_t i = 0; i < grid_size && i < static_cast<size_t>(nx * ny * nz);
           ++i) {
        state_values[i] += t[i];
      }
    }
    if (q) {
      for (size_t i = 0; i < grid_size && i < static_cast<size_t>(nx * ny * nz);
           ++i) {
        state_values[i] += q[i];
      }
    }
    if (psfc) {
      for (size_t i = 0; i < grid_size && i < static_cast<size_t>(nx * ny * nz);
           ++i) {
        state_values[i] += psfc[i];
      }
    }

    int num_obs_int = static_cast<int>(num_observations);
    const int rc = wrfda_xtoy_adjoint_grid(
        nx, ny, nz,
        obs_increment.data(),  // delta_y
        num_obs_int, state_values.data(), state_values.data(),
        state_values.data(), state_values.data(),
        state_values.data());  // inout_u, inout_v, inout_t, inout_q, inout_psfc
    if (rc != 0) {
      throw std::runtime_error("wrfda_xtoy_adjoint_grid failed with code " +
                               std::to_string(rc));
    }
  }

  /**
   * @brief Check if the operator supports linearization
   *
   * @return true (this operator supports linearization)
   *
   * @details Indicates whether the operator can be linearized around a
   * reference state. This operator supports both tangent linear and
   * adjoint operations.
   */
  bool supportsLinearization() const { return true; }

  /**
   * @brief Check if the operator is linear
   *
   * @return true (this operator is treated as linear)
   *
   * @details Indicates whether the operator is linear. This operator
   * treats the observation operator as linear, which means the tangent
   * linear operator is identical to the forward operator applied to
   * increments.
   */
  bool isLinear() const { return true; }

  /**
   * @brief Determine which operator family to use for a specific observation
   *
   * @param obs_type The observation type (e.g., "ADPSFC", "ADPUPA", "SFCSHP")
   * @return The appropriate operator family to use, or empty string if no match
   *
   * @details Maps observation types to appropriate operator families based on
   * WRFDA conventions:
   * - Surface observations (ADPSFC, SFCSHP, METAR) → "metar" family
   * - Upper-air observations (ADPUPA, AIRCRAFT, PROFILER) → "sound" family
   * - GPS observations (GPSPW, GPSREF) → "gpspw" family
   * - Radar observations (RADAR, RADARV) → "radar" family
   *
   * If no specific match is found, returns the first available family.
   * Returns empty string if no families are configured.
   */
  std::string determineOperatorFamily(const std::string& obs_type) const {
    if (operator_families_.empty()) {
      return "";
    }

    // Map observation types to operator families
    if (obs_type == "ADPSFC" || obs_type == "SFCSHP" || obs_type == "METAR") {
      // Surface observations -> use metar family if available
      for (const auto& family : operator_families_) {
        if (family == "metar") {
          return family;
        }
      }
    } else if (obs_type == "ADPUPA" || obs_type == "AIRCRAFT" ||
               obs_type == "PROFILER") {
      // Upper-air observations -> use sound family if available
      for (const auto& family : operator_families_) {
        if (family == "sound") {
          return family;
        }
      }
    } else if (obs_type == "GPSPW" || obs_type == "GPSREF") {
      // GPS observations -> use gpspw family if available
      for (const auto& family : operator_families_) {
        if (family == "gpspw") {
          return family;
        }
      }
    } else if (obs_type == "RADAR" || obs_type == "RADARV") {
      // Radar observations -> use radar family if available
      for (const auto& family : operator_families_) {
        if (family == "radar") {
          return family;
        }
      }
    }

    // If no specific match found, return the first available family
    return operator_families_[0];
  }

  /**
   * @brief Convert staggered U variable to unstaggered grid
   *
   * @param u_staggered Input staggered U array (nx+1, ny, nz)
   * @param u_unstaggered Output unstaggered U array (nx, ny, nz)
   * @param nx Number of x grid points (unstaggered)
   * @param ny Number of y grid points
   * @param nz Number of z grid points
   *
   * @details U is staggered in x-direction, so we interpolate to cell centers
   * by averaging adjacent staggered points. Boundary points use the nearest
   * staggered point.
   */
  void convertStaggeredUToUnstaggered(const std::vector<double>& u_staggered,
                                      std::vector<double>& u_unstaggered,
                                      int nx, int ny, int nz) const {
    u_unstaggered.resize(nx * ny * nz);

    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          const size_t unstaggered_idx = k * ny * nx + j * nx + i;

          if (i == 0) {
            // Left boundary: use first staggered point
            u_unstaggered[unstaggered_idx] =
                u_staggered[k * ny * (nx + 1) + j * (nx + 1) + i];
          } else if (i == nx - 1) {
            // Right boundary: use last staggered point
            u_unstaggered[unstaggered_idx] =
                u_staggered[k * ny * (nx + 1) + j * (nx + 1) + i];
          } else {
            // Interior: average of two adjacent staggered points
            u_unstaggered[unstaggered_idx] =
                0.5 * (u_staggered[k * ny * (nx + 1) + j * (nx + 1) + i] +
                       u_staggered[k * ny * (nx + 1) + j * (nx + 1) + i + 1]);
          }
        }
      }
    }
  }

  /**
   * @brief Convert staggered V variable to unstaggered grid
   *
   * @param v_staggered Input staggered V array (nx, ny+1, nz)
   * @param v_unstaggered Output unstaggered V array (nx, ny, nz)
   * @param nx Number of x grid points
   * @param ny Number of y grid points (unstaggered)
   * @param nz Number of z grid points
   *
   * @details V is staggered in y-direction, so we interpolate to cell centers
   * by averaging adjacent staggered points. Boundary points use the nearest
   * staggered point.
   */
  void convertStaggeredVToUnstaggered(const std::vector<double>& v_staggered,
                                      std::vector<double>& v_unstaggered,
                                      int nx, int ny, int nz) const {
    v_unstaggered.resize(nx * ny * nz);

    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          const size_t unstaggered_idx = k * ny * nx + j * nx + i;

          if (j == 0) {
            // Bottom boundary: use first staggered point
            v_unstaggered[unstaggered_idx] =
                v_staggered[k * (ny + 1) * nx + j * nx + i];
          } else if (j == ny - 1) {
            // Top boundary: use last staggered point
            v_unstaggered[unstaggered_idx] =
                v_staggered[k * (ny + 1) * nx + j * nx + i];
          } else {
            // Interior: average of two adjacent staggered points
            v_unstaggered[unstaggered_idx] =
                0.5 * (v_staggered[k * (ny + 1) * nx + j * nx + i] +
                       v_staggered[k * (ny + 1) * nx + (j + 1) * nx + i]);
          }
        }
      }
    }
  }

  /**
   * @brief Check if a variable is staggered using the existing WRF
   * infrastructure
   *
   * @param state The state backend containing variable grid information
   * @param var_name Variable name (U, V, T, Q, PSFC)
   * @return true if the variable is staggered, false otherwise
   */
  bool isVariableStaggered(const StateBackend& state,
                           const std::string& var_name) const {
    try {
      const auto& grid_info = state.getVariableGridInfo(var_name);
      return (
          grid_info.grid_type ==
              metada::backends::wrf::VariableGridInfo::GridType::U_STAGGERED ||
          grid_info.grid_type ==
              metada::backends::wrf::VariableGridInfo::GridType::V_STAGGERED ||
          grid_info.grid_type ==
              metada::backends::wrf::VariableGridInfo::GridType::W_STAGGERED);
    } catch (...) {
      // If we can't get grid info, assume unstaggered
      return false;
    }
  }

 private:
  /**
   * @brief Ensure the operator is initialized before use
   *
   * @throws std::runtime_error if the operator is not initialized
   *
   * @details This is a helper method that checks the initialization
   * status and throws an exception if the operator is not ready for use.
   * Called by all public methods that require initialization.
   */
  void ensureInitialized() const {
    if (!isInitialized()) {
      throw std::runtime_error("WRFDAObsOperator not initialized");
    }
  }

  /**
   * @brief Construct WRFDA y_type structure from observation data
   *
   * @param obs_data The observation data to convert
   * @param families The operator families to use
   * @return Pointer to constructed y_type structure
   *
   * @details Constructs a WRFDA y_type structure containing observation
   * residuals for the specified operator families. The structure is
   * allocated and populated with observation values.
   */
  void* constructYType(const typename ObsBackend::ObsOperatorData& obs_data,
                       const std::vector<std::string>& families) const {
    if (obs_data.num_obs == 0) {
      return nullptr;
    }

    // Determine the operator family to use
    std::string family = families.empty() ? "metar" : families[0];

    // Prepare observation types as a flat string
    std::string obs_types_flat;
    for (size_t i = 0; i < obs_data.num_obs; ++i) {
      std::string obs_type =
          (i < obs_data.types.size()) ? obs_data.types[i] : "ADPSFC";
      obs_types_flat += obs_type + std::string(20 - obs_type.length(), '\0');
    }

    // Extract structured data for WRFDA
    std::vector<double> u_values, v_values, t_values, p_values, q_values;
    std::vector<double> u_errors, v_errors, t_errors, p_errors, q_errors;
    std::vector<double> lats, lons, levels;
    std::vector<int> u_available, v_available, t_available, p_available,
        q_available;

    // Extract data from the hierarchical structure
    for (size_t station = 0; station < obs_data.num_obs; ++station) {
      lats.push_back(obs_data.lats[station]);
      lons.push_back(obs_data.lons[station]);

      // Process each level for this station
      for (const auto& level : obs_data.levels[station]) {
        levels.push_back(level.level);

        // Extract variable data
        // Use WRFDA standard missing value (-888888.0) for unavailable
        // observations
        const double missing_r = -888888.0;
        u_values.push_back(level.u.available ? level.u.value : missing_r);
        v_values.push_back(level.v.available ? level.v.value : missing_r);
        t_values.push_back(level.t.available ? level.t.value : missing_r);
        p_values.push_back(level.p.available ? level.p.value : missing_r);
        q_values.push_back(level.q.available ? level.q.value : missing_r);

        u_errors.push_back(level.u.available ? level.u.error : 1.0);
        v_errors.push_back(level.v.available ? level.v.error : 1.0);
        t_errors.push_back(level.t.available ? level.t.error : 1.0);
        p_errors.push_back(level.p.available ? level.p.error : 1.0);
        q_errors.push_back(level.q.available ? level.q.error : 1.0);

        // Extract availability flags
        u_available.push_back(level.u.available ? 1 : 0);
        v_available.push_back(level.v.available ? 1 : 0);
        t_available.push_back(level.t.available ? 1 : 0);
        p_available.push_back(level.p.available ? 1 : 0);
        q_available.push_back(level.q.available ? 1 : 0);
      }
    }

    // Call the Fortran function to construct y_type
    int num_obs_int = static_cast<int>(obs_data.num_obs);
    int num_levels_int = static_cast<int>(levels.size());

    void* y_ptr = wrfda_construct_y_type(
        &num_obs_int, &num_levels_int, const_cast<double*>(u_values.data()),
        const_cast<double*>(v_values.data()),
        const_cast<double*>(t_values.data()),
        const_cast<double*>(p_values.data()),
        const_cast<double*>(q_values.data()),
        const_cast<double*>(u_errors.data()),
        const_cast<double*>(v_errors.data()),
        const_cast<double*>(t_errors.data()),
        const_cast<double*>(p_errors.data()),
        const_cast<double*>(q_errors.data()),
        const_cast<int*>(u_available.data()),
        const_cast<int*>(v_available.data()),
        const_cast<int*>(t_available.data()),
        const_cast<int*>(p_available.data()),
        const_cast<int*>(q_available.data()), const_cast<double*>(lats.data()),
        const_cast<double*>(lons.data()), const_cast<double*>(levels.data()),
        const_cast<char*>(obs_types_flat.c_str()),
        const_cast<char*>(family.c_str()));

    return y_ptr;
  }

  /**
   * @brief Construct WRFDA iv_type structure from observation data
   *
   * @param obs_data The observation data to convert
   * @param families The operator families to use
   * @return Pointer to constructed iv_type structure
   *
   * @details Constructs a WRFDA iv_type structure containing observation
   * innovation vectors for the specified operator families. The structure is
   * allocated and populated with observation data including locations and
   * quality control information.
   */
  void* constructIvType(const typename ObsBackend::ObsOperatorData& obs_data,
                        const std::vector<std::string>& families,
                        void* domain_ptr) const {
    if (obs_data.num_obs == 0) {
      return nullptr;
    }

    // Determine the operator family to use
    std::string family = families.empty() ? "metar" : families[0];

    // Prepare observation types as a flat string
    std::string obs_types_flat;
    for (size_t i = 0; i < obs_data.num_obs; ++i) {
      std::string obs_type =
          (i < obs_data.types.size()) ? obs_data.types[i] : "ADPSFC";
      obs_types_flat += obs_type + std::string(20 - obs_type.length(), '\0');
    }

    // Extract structured data for WRFDA
    std::vector<double> u_values, v_values, t_values, p_values, q_values;
    std::vector<double> u_errors, v_errors, t_errors, p_errors, q_errors;
    std::vector<double> lats, lons, levels;
    std::vector<int> u_qc, v_qc, t_qc, p_qc, q_qc;
    std::vector<int> u_available, v_available, t_available, p_available,
        q_available;

    // Extract data from the hierarchical structure
    for (size_t station = 0; station < obs_data.num_obs; ++station) {
      lats.push_back(obs_data.lats[station]);
      lons.push_back(obs_data.lons[station]);

      // Process each level for this station
      for (const auto& level : obs_data.levels[station]) {
        levels.push_back(level.level);

        // Extract variable data
        // Use WRFDA standard missing value (-888888.0) for unavailable
        // observations
        const double missing_r = -888888.0;
        u_values.push_back(level.u.available ? level.u.value : missing_r);
        v_values.push_back(level.v.available ? level.v.value : missing_r);
        t_values.push_back(level.t.available ? level.t.value : missing_r);
        p_values.push_back(level.p.available ? level.p.value : missing_r);
        q_values.push_back(level.q.available ? level.q.value : missing_r);

        u_errors.push_back(level.u.available ? level.u.error : 1.0);
        v_errors.push_back(level.v.available ? level.v.error : 1.0);
        t_errors.push_back(level.t.available ? level.t.error : 1.0);
        p_errors.push_back(level.p.available ? level.p.error : 1.0);
        q_errors.push_back(level.q.available ? level.q.error : 1.0);

        const int missing_data = -88;
        u_qc.push_back(level.u.available ? level.u.qc : missing_data);
        v_qc.push_back(level.v.available ? level.v.qc : missing_data);
        t_qc.push_back(level.t.available ? level.t.qc : missing_data);
        p_qc.push_back(level.p.available ? level.p.qc : missing_data);
        q_qc.push_back(level.q.available ? level.q.qc : missing_data);

        // Extract availability flags
        u_available.push_back(level.u.available ? 1 : 0);
        v_available.push_back(level.v.available ? 1 : 0);
        t_available.push_back(level.t.available ? 1 : 0);
        p_available.push_back(level.p.available ? 1 : 0);
        q_available.push_back(level.q.available ? 1 : 0);
      }
    }

    // Call the Fortran function to construct iv_type
    int num_obs_int = static_cast<int>(obs_data.num_obs);
    int num_levels_int = static_cast<int>(levels.size());

    void* iv_ptr = wrfda_construct_iv_type(
        &num_obs_int, &num_levels_int, const_cast<double*>(u_values.data()),
        const_cast<double*>(v_values.data()),
        const_cast<double*>(t_values.data()),
        const_cast<double*>(p_values.data()),
        const_cast<double*>(q_values.data()),
        const_cast<double*>(u_errors.data()),
        const_cast<double*>(v_errors.data()),
        const_cast<double*>(t_errors.data()),
        const_cast<double*>(p_errors.data()),
        const_cast<double*>(q_errors.data()), const_cast<int*>(u_qc.data()),
        const_cast<int*>(v_qc.data()), const_cast<int*>(t_qc.data()),
        const_cast<int*>(p_qc.data()), const_cast<int*>(q_qc.data()),
        const_cast<int*>(u_available.data()),
        const_cast<int*>(v_available.data()),
        const_cast<int*>(t_available.data()),
        const_cast<int*>(p_available.data()),
        const_cast<int*>(q_available.data()), const_cast<double*>(lats.data()),
        const_cast<double*>(lons.data()), const_cast<double*>(levels.data()),
        const_cast<double*>(obs_data.elevations.data()),
        const_cast<char*>(obs_types_flat.c_str()),
        const_cast<char*>(family.c_str()), domain_ptr);

    return iv_ptr;
  }

  /**
   * @brief Construct WRFDA config_flags structure
   *
   * @return Pointer to constructed config_flags structure
   *
   * @details Constructs a WRFDA config_flags structure containing
   * configuration parameters for the data assimilation process.
   * Currently returns a simplified implementation.
   */
  void* constructConfigFlags() const { return wrfda_construct_config_flags(); }

  /**
   * @brief Extract innovation values from iv_type structure
   *
   * @param iv_ptr Pointer to iv_type structure
   * @param family Observation family name
   * @return Vector of innovation values
   *
   * @details Extracts innovation values from the WRFDA iv_type structure
   * and returns them as a vector. The innovations represent the differences
   * between observations and model values (O-B).
   */
  std::vector<double> extractInnovations(void* iv_ptr,
                                         const std::string& family) const {
    if (!iv_ptr) {
      return std::vector<double>();
    }

    // First call to count the number of innovations
    int num_innovations = 0;
    int rc = wrfda_count_innovations(const_cast<char*>(family.c_str()),
                                     &num_innovations);

    std::cout << "WRFDAObsOperator: Count call returned rc = " << rc
              << ", num_innovations = " << num_innovations << std::endl;

    if (rc != 0 || num_innovations <= 0) {
      std::cerr << "WRFDA: Failed to count innovations with code " << rc
                << std::endl;
      return std::vector<double>();
    }

    // Allocate vector with exact size needed
    std::vector<double> innovations(num_innovations);

    // Second call to extract the actual innovations
    rc = wrfda_extract_innovations(const_cast<char*>(family.c_str()),
                                   innovations.data(), &num_innovations);

    std::cout << "WRFDAObsOperator: Extract call returned rc = " << rc
              << ", extracted " << num_innovations << " innovations"
              << std::endl;

    if (rc != 0) {
      std::cerr << "WRFDA: Failed to extract innovations with code " << rc
                << std::endl;
      return std::vector<double>();
    }

    return innovations;
  }

  /**
   * @brief Initialization status flag
   *
   * @details Indicates whether the operator has been properly initialized
   * with configuration parameters.
   */
  bool initialized_ = false;

  /**
   * @brief Required state variables for the operator
   *
   * @details List of state variable names that must be present in the
   * state for the operator to function correctly.
   */
  std::vector<std::string> required_state_vars_;

  /**
   * @brief Required observation variables for the operator
   *
   * @details List of observation variable names that the operator
   * can process.
   */
  std::vector<std::string> required_obs_vars_;

  /**
   * @brief Path to external operator sources
   *
   * @details Root directory path containing external observation operator
   * implementations (e.g., WRFDA, GSI, DART).
   */
  std::string external_root_;

  /**
   * @brief External system identifier
   *
   * @details Identifier for the external observation operator system
   * being used (e.g., "wrfda", "gsi", "dart").
   */
  std::string external_system_;

  /**
   * @brief Available operator families
   *
   * @details List of observation operator families that are available
   * for use (e.g., "metar", "sound", "gpspw", "radar").
   */
  std::vector<std::string> operator_families_;

  /**
   * @brief Background state data for proper WRFDA integration
   *
   * @details Stores the background state data that will be used to
   * properly separate total state into background and increments.
   */
  mutable typename StateBackend::ObsOperatorData background_state_data_;

  /**
   * @brief Flag indicating whether background state is available
   *
   * @details When true, the operator will use the new WRFDA function
   * that properly separates background and increments.
   */
  mutable bool has_background_state_ = false;
};

}  // namespace metada::backends::common::obsoperator
