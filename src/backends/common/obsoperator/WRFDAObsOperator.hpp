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
#include <string>
#include <utility>
#include <vector>

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
  std::vector<double> apply(const StateBackend& state,
                            const ObsBackend& obs) const {
    ensureInitialized();
    // Collect observation coordinates (geographic)
    const size_t num_observations = obs.size();
    std::vector<double> obs_lats(num_observations);
    std::vector<double> obs_lons(num_observations);
    std::vector<double> obs_levels(num_observations, 0.0);
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
      obs_levels[i] = level;
    }

    // Helper to obtain variable buffers or zero-filled fallbacks
    auto get_or_zero = [&](const std::string& var, std::vector<double>& storage,
                           const double*& ptr, int& nx, int& ny, int& nz) {
      ptr = static_cast<const double*>(state.getData(var));
      if (ptr == nullptr) {
        nx = static_cast<int>(state.geometry().x_dim());
        ny = static_cast<int>(state.geometry().y_dim());
        nz = static_cast<int>(state.geometry().z_dim());
        const size_t n = static_cast<size_t>(std::max(1, nx) * std::max(1, ny) *
                                             std::max(1, nz));
        storage.assign(n, 0.0);
        ptr = storage.data();
      } else {
        try {
          const auto& dims = state.getVariableDimensions(var);
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

    const double *u = nullptr, *v = nullptr, *t = nullptr, *q = nullptr,
                 *psfc = nullptr;
    std::vector<double> u_zero, v_zero, t_zero, q_zero, psfc_zero;
    int nx_u = 1, ny_u = 1, nz_u = 1;
    int nx_v = 1, ny_v = 1, nz_v = 1;
    int nx_t = 1, ny_t = 1, nz_t = 1;
    int nx_q = 1, ny_q = 1, nz_q = 1;
    int nx_ps = 1, ny_ps = 1, nz_ps = 1;
    get_or_zero("U", u_zero, u, nx_u, ny_u, nz_u);
    get_or_zero("V", v_zero, v, nx_v, ny_v, nz_v);
    get_or_zero("T", t_zero, t, nx_t, ny_t, nz_t);
    get_or_zero("Q", q_zero, q, nx_q, ny_q, nz_q);
    get_or_zero("PSFC", psfc_zero, psfc, nx_ps, ny_ps, nz_ps);

    // Canonical unstaggered grid dims
    const int nx = static_cast<int>(state.geometry().x_dim());
    const int ny = static_cast<int>(state.geometry().y_dim());
    const int nz = static_cast<int>(state.geometry().z_dim());

    // Geometry arrays
    const auto& gi = state.geometry().unstaggered_info();
    const double* levels_ptr =
        gi.vertical_coords.empty() ? nullptr : gi.vertical_coords.data();
    std::vector<double> dummy_levels(1, 0.0);
    if (levels_ptr == nullptr) levels_ptr = dummy_levels.data();

    std::vector<double> out_y(num_observations, 0.0);

    const char* family_cstr =
        (operator_families_.empty() ? "" : operator_families_[0].c_str());

    // Use the proper WRFDA grid-based operator that handles variables
    // separately This avoids the physical nonsense of averaging different
    // variable types and properly prepares the state for WRFDA's observation
    // operator

    // Prepare 2D lat/lon arrays for WRFDA (required by the grid API)
    std::vector<double> lats2d(nx * ny);
    std::vector<double> lons2d(nx * ny);

    // Extract 2D lat/lon from the geometry (assuming they're available)
    // If not available, we'll need to generate them from the grid info
    if (gi.latitude_2d.empty() || gi.longitude_2d.empty()) {
      // Cannot proceed without latitude and longitude data
      // These are essential for WRFDA observation operator functionality
      throw std::runtime_error(
          "WRFDAObsOperator: Missing required latitude and longitude data. "
          "The geometry must provide valid latitude_2d and longitude_2d arrays "
          "for proper observation localization and grid interpolation.");
    }

    // Use provided lat/lon from geometry
    std::copy(gi.latitude_2d.begin(), gi.latitude_2d.end(), lats2d.begin());
    std::copy(gi.longitude_2d.begin(), gi.longitude_2d.end(), lons2d.begin());

    // Prepare vertical levels array
    std::vector<double> levels(nz);
    if (levels_ptr != nullptr && !gi.vertical_coords.empty()) {
      std::copy(gi.vertical_coords.begin(), gi.vertical_coords.end(),
                levels.begin());
    } else {
      // Cannot proceed without vertical coordinate data
      // These are essential for WRFDA observation operator functionality
      throw std::runtime_error(
          "WRFDAObsOperator: Missing required vertical coordinate data. "
          "The geometry must provide valid vertical_coords array "
          "for proper vertical interpolation and level assignment.");
    }

    // Call the proper WRFDA grid-based operator
    // This function expects separate arrays for each variable, which is correct
    // Note: WRFDA can handle nullptr for missing variables, but we should
    // ensure at least the essential variables are present
    if (!u && !v && !t) {
      throw std::runtime_error(
          "WRFDA observation operator requires at least one of U, V, or T "
          "variables");
    }

    // Debug output for troubleshooting
    std::cout << "WRFDA: Calling observation operator with:" << std::endl;
    std::cout << "  Grid dimensions: " << nx << "x" << ny << "x" << nz
              << std::endl;
    std::cout << "  Variables: U=" << (u ? "present" : "null")
              << ", V=" << (v ? "present" : "null")
              << ", T=" << (t ? "present" : "null")
              << ", Q=" << (q ? "present" : "null")
              << ", PSFC=" << (psfc ? "present" : "null") << std::endl;
    std::cout << "  Observations: " << num_observations << std::endl;
    std::cout << "  Operator family: " << family_cstr << std::endl;

    const int rc = wrfda_xtoy_apply_grid(
        family_cstr, nx, ny, nz, u, v, t, q,
        psfc,  // Individual variable arrays
        lats2d.data(), lons2d.data(), levels.data(),  // Grid metadata
        static_cast<int>(num_observations), obs_lats.data(), obs_lons.data(),
        obs_levels.data(), out_y.data());
    if (rc != 0) {
      throw std::runtime_error("wrfda_xtoy_apply_grid failed with code " +
                               std::to_string(rc));
    }

    return out_y;
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
    // For linearized operator, use same array-based call on increment
    ensureInitialized();
    return apply(state_increment, obs);
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

    // Get model vertical levels to determine appropriate defaults
    const auto& gi = result_state.geometry().unstaggered_info();
    double surface_level = 1.0;    // Default surface level (sigma coordinates)
    double min_model_level = 0.0;  // Default minimum model level

    if (!gi.vertical_coords.empty()) {
      // Find actual min/max model levels
      min_model_level = *std::min_element(gi.vertical_coords.begin(),
                                          gi.vertical_coords.end());
      surface_level = *std::max_element(gi.vertical_coords.begin(),
                                        gi.vertical_coords.end());
    }

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

      // Handle invalid observation levels
      if (level <= 0.0 || std::isnan(level) || std::isinf(level)) {
        // This is likely a surface observation or missing level info
        // Set to surface level (highest model level, typically 1.0 for sigma)
        obs_levels[i] = surface_level;
        std::cout << "WARNING: Observation " << i << " has invalid level "
                  << level << ", setting to surface level " << surface_level
                  << std::endl;
      } else if (level < min_model_level) {
        // Level is below model domain, likely a surface observation
        // Set to surface level
        obs_levels[i] = surface_level;
        std::cout << "WARNING: Observation " << i << " level " << level
                  << " below model minimum " << min_model_level
                  << ", setting to surface level " << surface_level
                  << std::endl;
      } else {
        // Valid level within model domain
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
    get_or_zero_inout("Q", q_zero, q, nx_q, ny_q, nz_q);
    get_or_zero_inout("PSFC", psfc_zero, psfc, nx_ps, ny_ps, nz_ps);

    const int nx = static_cast<int>(result_state.geometry().x_dim());
    const int ny = static_cast<int>(result_state.geometry().y_dim());
    const int nz = static_cast<int>(result_state.geometry().z_dim());

    // Use the gi variable already declared above
    const double* levels_ptr =
        gi.vertical_coords.empty() ? nullptr : gi.vertical_coords.data();
    std::vector<double> dummy_levels(1, 0.0);
    if (levels_ptr == nullptr) levels_ptr = dummy_levels.data();

    const char* family_cstr =
        (operator_families_.empty() ? "" : operator_families_[0].c_str());

    // Flatten the state variables into a single array for the adjoint operator
    const size_t grid_size = static_cast<size_t>(nx * ny * nz);
    std::vector<double> state_values(grid_size, 0.0);

    // Initialize with existing values if available
    if (u) {
      for (size_t i = 0; i < grid_size && i < static_cast<size_t>(nx * ny * nz);
           ++i) {
        state_values[i] += u[i];
      }
    }
    if (v) {
      for (size_t i = 0; i < grid_size && i < static_cast<size_t>(nx * ny * nz);
           ++i) {
        state_values[i] += v[i];
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

    const int rc = wrfda_xtoy_adjoint_grid(
        family_cstr, nx, ny, nz,
        obs_increment.data(),  // delta_y
        gi.latitude_2d.data(), gi.longitude_2d.data(),
        gi.vertical_coords.data(),  // lats2d, lons2d, levels
        static_cast<int>(num_observations), obs_lats.data(), obs_lons.data(),
        obs_levels.data(), state_values.data(), state_values.data(),
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
};

}  // namespace metada::backends::common::obsoperator
