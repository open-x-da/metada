/**
 * @file MACOMObsOperator.hpp
 * @brief MACOM observation operator backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This class provides an observation operator specifically designed for MACOM's
 * unstructured hexagonal grid. It uses the existing getValuesAtNearestPoints
 * method for efficient spatial interpolation between geographic observation
 * locations and the irregular MACOM grid points.
 */

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "Location.hpp"
#include "PointObservation.hpp"

namespace metada::backends::macom {

using framework::CoordinateSystem;
using framework::Location;
using framework::ObservationPoint;

/**
 * @brief MACOM observation operator backend implementation
 *
 * @details
 * This class implements an observation operator specifically for MACOM that:
 * - Maps state variables from unstructured hexagonal grid to observation
 * locations
 * - Uses existing getValuesAtNearestPoints for efficient spatial searching
 * - Supports geographic coordinate observations
 * - Works with MACOM's native grid structure
 *
 * @tparam StateBackend Type of MACOM state backend
 * @tparam ObsBackend Type of observation backend
 */
template <typename StateBackend, typename ObsBackend>
class MACOMObsOperator {
 public:
  // =============================================================================
  // FRAMEWORK CONCEPTS REQUIRED INTERFACES
  // Required by ObsOperatorBackendImpl concept
  // =============================================================================

  // Type aliases for concept compliance
  using iterator_type = typename std::vector<double>::iterator;
  using const_iterator_type = typename std::vector<double>::const_iterator;

  // --- Resource management (required by framework) ---

  /**
   * @brief Default constructor is deleted (required by framework)
   */
  MACOMObsOperator() = delete;

  /**
   * @brief Copy constructor is deleted (required by framework)
   */
  MACOMObsOperator(const MACOMObsOperator&) = delete;

  /**
   * @brief Copy assignment operator is deleted (required by framework)
   */
  MACOMObsOperator& operator=(const MACOMObsOperator&) = delete;

  /**
   * @brief Constructor that initializes from configuration (required by
   * framework)
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  explicit MACOMObsOperator(const ConfigBackend& config)
      : interpolation_method_("nearest_neighbor") {
    // Read configuration parameters
    try {
      // Optionally read interpolation method
      if (config.HasKey("interpolation_method")) {
        interpolation_method_ = config.Get("interpolation_method").asString();
      }

      // Read required variables
      const auto& variables_configs = config.Get("variables").asVectorMap();
      for (const auto& var_map : variables_configs) {
        for (const auto& [var_name, var_config] : var_map) {
          const auto& var_backend = ConfigBackend(var_config.asMap());
          bool var_if_use = var_backend.Get("if_use").asBool();
          if (var_if_use) {
            required_state_vars_.push_back(var_name);
            required_obs_vars_.push_back(var_name);
          }
        }
      }

      initialized_ = true;
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to initialize MACOMObsOperator: " +
                               std::string(e.what()));
    }
  }

  /**
   * @brief Move constructor (required by framework)
   */
  MACOMObsOperator(MACOMObsOperator&& other) noexcept = default;

  /**
   * @brief Move assignment operator (required by framework)
   */
  MACOMObsOperator& operator=(MACOMObsOperator&& other) noexcept = default;

  /**
   * @brief Destructor (required by framework)
   */
  ~MACOMObsOperator() = default;

  /**
   * @brief Clone this observation operator (required by framework)
   * @return A new instance that is a copy of this one
   */
  MACOMObsOperator clone() const {
    MACOMObsOperator new_op(*this, true);
    return new_op;
  }

  // --- Initialization interface (required by framework) ---

  /**
   * @brief Check if the observation operator is initialized
   * @return true if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Initialize the observation operator (required by framework)
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  void initialize(const ConfigBackend& config) {
    try {
      // Optionally read interpolation method
      if (config.HasKey("interpolation_method")) {
        interpolation_method_ = config.Get("interpolation_method").asString();
      }

      // Read required variables
      if (config.HasKey("variables")) {
        const auto& variables_configs = config.Get("variables").asVectorMap();
        for (const auto& var_map : variables_configs) {
          for (const auto& [var_name, var_config] : var_map) {
            const auto& var_backend = ConfigBackend(var_config.asMap());
            bool var_if_use = var_backend.Get("if_use").asBool();
            if (var_if_use) {
              required_state_vars_.push_back(var_name);
              required_obs_vars_.push_back(var_name);
            }
          }
        }
      }

      initialized_ = true;
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to initialize MACOMObsOperator: " +
                               std::string(e.what()));
    }
  }

  /**
   * @brief Initialize the observation operator (no-param version)
   */
  void initialize() {
    if (!initialized_) {
      initialized_ = true;
    }
  }

  // --- Core operator interface (required by framework) ---

  /**
   * @brief Apply observation operator: H(x) -> y_obs (required by framework)
   *
   * Maps state variables to observation space using MACOM's unstructured grid.
   * This uses the existing getValuesAtNearestPoints method for efficient
   * spatial interpolation.
   *
   * @param state Model state on MACOM grid
   * @param obs Observations with geographic locations
   * @return Vector of interpolated values at observation locations
   */
  std::vector<double> apply(const StateBackend& state,
                            const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("MACOMObsOperator not initialized");
    }

    std::vector<double> result;
    result.reserve(obs.size());

    // Prepare observation coordinates for batch processing
    std::vector<double> obs_lats, obs_lons, obs_levels;
    obs_lats.reserve(obs.size());
    obs_lons.reserve(obs.size());
    obs_levels.reserve(obs.size());

    for (const auto& obs_point : obs) {
      if (!obs_point.is_valid) {
        obs_lats.push_back(0.0);  // Dummy values for invalid observations
        obs_lons.push_back(0.0);
        obs_levels.push_back(0.0);
      } else if (obs_point.location.getCoordinateSystem() ==
                 CoordinateSystem::GEOGRAPHIC) {
        auto [lat, lon, level] = obs_point.location.getGeographicCoords();
        obs_lats.push_back(lat);
        obs_lons.push_back(lon);
        obs_levels.push_back(level);
      } else {
        throw std::runtime_error(
            "MACOMObsOperator only supports geographic coordinates for "
            "observations");
      }
    }

    // Use the existing getValuesAtNearestPoints method for efficient batch
    // processing Use active variable and enable horizontal interpolation for
    // better accuracy
    std::string active_var = state.getActiveVariable();
    auto interpolated_values = state.getValuesAtNearestPoints(
        obs_lons, obs_lats, obs_levels, active_var, false);

    // Process results and handle invalid observations
    size_t obs_idx = 0;
    for (const auto& obs_point : obs) {
      if (!obs_point.is_valid) {
        result.push_back(0.0);  // Invalid observations get zero
      } else {
        result.push_back(interpolated_values[obs_idx]);
      }
      obs_idx++;
    }

    return result;
  }

  /**
   * @brief Apply adjoint observation operator (required by framework)
   *
   * Maps observation increments back to state space.
   *
   * @param obs_increment Observation increment
   * @param reference_state Reference state
   * @param increment_state Increment state (output)
   * @param obs Observations
   */
  void applyAdjoint(const std::vector<double>& obs_increment,
                    const StateBackend& reference_state,
                    StateBackend& increment_state,
                    const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("MACOMObsOperator not initialized");
    }

    // For linear observation operator, adjoint is the transpose
    // This is a simplified implementation - can be extended for more complex
    // cases
    // Use move assignment to copy the structure from reference state
    increment_state = std::move(*(reference_state.clone()));
    increment_state.zero();

    // Simple adjoint implementation: distribute observation increments back to
    // grid This is a placeholder - should be implemented based on the actual
    // adjoint
    for (size_t i = 0; i < obs_increment.size() && i < obs.size(); ++i) {
      const auto& obs_point = obs[i];
      if (obs_point.is_valid) {
        // For now, just add the increment to the nearest grid point
        // This should be replaced with proper adjoint interpolation
        if (obs_point.location.getCoordinateSystem() ==
            CoordinateSystem::GEOGRAPHIC) {
          auto [lat, lon, level] = obs_point.location.getGeographicCoords();
          increment_state.modifyValueAtLocation(lat, lon, obs_increment[i]);
        }
      }
    }
  }

  /**
   * @brief Apply tangent linear observation operator: H dx
   *
   * @details For MACOM, the tangent linear operator is identical to the forward
   * operator since we use linear interpolation (nearest neighbor).
   *
   * @param state_increment State increment to transform
   * @param reference_state Reference state around which to linearize
   * @param obs Reference observations for context
   * @return Vector containing the transformed increment in observation space
   */
  std::vector<double> applyTangentLinear(const StateBackend& state_increment,
                                         const StateBackend& reference_state,
                                         const ObsBackend& obs) const {
    // For linear interpolation, tangent linear is identical to forward operator
    return apply(state_increment, obs);
  }

  /**
   * @brief Check if tangent linear and adjoint operators are available
   *
   * @return True if tangent linear/adjoint operators are supported
   */
  bool supportsLinearization() const {
    return true;  // MACOM supports both TL and AD
  }

  /**
   * @brief Check if the observation operator is linear
   *
   * @details MACOM uses linear interpolation (nearest neighbor), so the
   * operator is linear.
   *
   * @return True if the observation operator is linear
   */
  bool isLinear() const {
    return true;  // MACOM uses linear interpolation
  }

  // --- Variable interface (required by framework) ---

  /**
   * @brief Get the state variables required by this observation operator
   * @return Const reference to vector of required state variable names
   */
  const std::vector<std::string>& getRequiredStateVars() const {
    return required_state_vars_;
  }

  /**
   * @brief Get the observation variables required by this observation operator
   * @return Const reference to vector of required observation variable names
   */
  const std::vector<std::string>& getRequiredObsVars() const {
    return required_obs_vars_;
  }

  // =============================================================================
  // MACOM SPECIFIC FUNCTIONALITY
  // These are MACOM-specific methods beyond framework requirements
  // =============================================================================

  /**
   * @brief Get the interpolation method being used
   * @return String describing the interpolation method
   */
  const std::string& getInterpolationMethod() const {
    return interpolation_method_;
  }

  /**
   * @brief Set the interpolation method
   * @param method The interpolation method to use
   */
  void setInterpolationMethod(const std::string& method) {
    interpolation_method_ = method;
  }

 private:
  /**
   * @brief Private constructor for cloning
   */
  MACOMObsOperator(const MACOMObsOperator& other, bool)
      : interpolation_method_(other.interpolation_method_),
        required_state_vars_(other.required_state_vars_),
        required_obs_vars_(other.required_obs_vars_),
        initialized_(other.initialized_) {}

  std::string interpolation_method_;              ///< Interpolation method
  std::vector<std::string> required_state_vars_;  ///< Required state variables
  std::vector<std::string>
      required_obs_vars_;     ///< Required observation variables
  bool initialized_ = false;  ///< Initialization status
};

}  // namespace metada::backends::macom