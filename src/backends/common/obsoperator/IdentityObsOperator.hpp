/**
 * @file IdentityObsOperator.hpp
 * @brief Identity observation operator backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This class provides an identity observation operator backend interface that maps 
 * any state backend to GridObservation. The observation operator performs nearest-neighbor 
 * interpolation from the model grid to observation locations.
 */

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "GridObservation.hpp"
#include "Location.hpp"
#include "PointObservation.hpp"

namespace metada::backends::common::obsoperator {

using framework::CoordinateSystem;
using common::observation::GridObservation;

/**
 * @brief Identity observation operator backend implementation
 *
 * @details
 * This class implements an identity observation operator that:
 * - Maps state variables from model grid to observation locations
 * - Uses nearest-neighbor interpolation for grid-to-point mapping
 * - Supports multiple observation types and variables
 * - Handles missing values and quality control
 * - Provides configuration-based initialization
 * - Works with any state backend that provides required interface
 *
 * @tparam StateBackend Type of state backend to operate on
 */
template <typename StateBackend>
class IdentityObsOperator {
 public:
  // Delete default constructor
  IdentityObsOperator() = delete;

  // Delete copy constructor and assignment
  IdentityObsOperator(const IdentityObsOperator&) = delete;
  IdentityObsOperator& operator=(const IdentityObsOperator&) = delete;

  /**
   * @brief Constructor that initializes from configuration
   * @tparam ConfigBackend Type of configuration backend
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  explicit IdentityObsOperator(const ConfigBackend& config) {
    initialize(config);
  }

  /**
   * @brief Move constructor
   * @param other Observation operator to move from
   */
  IdentityObsOperator(IdentityObsOperator&& other) noexcept
      : initialized_(other.initialized_),
        required_state_vars_(std::move(other.required_state_vars_)),
        required_obs_vars_(std::move(other.required_obs_vars_)) {
    other.initialized_ = false;
  }

  /**
   * @brief Move assignment operator
   * @param other Observation operator to move from
   * @return Reference to this observation operator
   */
  IdentityObsOperator& operator=(IdentityObsOperator&& other) noexcept {
    if (this != &other) {
      initialized_ = other.initialized_;
      required_state_vars_ = std::move(other.required_state_vars_);
      required_obs_vars_ = std::move(other.required_obs_vars_);
      other.initialized_ = false;
    }
    return *this;
  }

  /**
   * @brief Initialize with configuration
   * @tparam ConfigBackend Type of configuration backend
   * @param config Configuration backend instance
   */
  template <typename ConfigBackend>
  void initialize(const ConfigBackend& config) {
    if (isInitialized()) {
      throw std::runtime_error("IdentityObsOperator already initialized");
    }

    // Get required state variables from config
    try {
      required_state_vars_ = config.Get("required_state_vars").asVectorString();
    } catch (...) {
      // Default to "state" if not specified
      required_state_vars_ = {"state"};
    }

    // Get required observation variables from config
    try {
      required_obs_vars_ = config.Get("required_obs_vars").asVectorString();
    } catch (...) {
      // Default to empty if not specified
      required_obs_vars_ = {};
    }

    initialized_ = true;
  }

  /**
   * @brief Check if the observation operator is initialized
   * @return True if initialized, false otherwise
   */
  bool isInitialized() const { return initialized_; }

  /**
   * @brief Apply forward observation operator: H(x)
   *
   * @details Maps state to observation space by interpolating from model grid
   * to observation locations using nearest-neighbor interpolation.
   *
   * @param state Model state to transform
   * @param obs Output observation to store the result
   * @return Vector of interpolated values at observation locations
   */
  std::vector<double> apply(const StateBackend& state,
                            const GridObservation& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("IdentityObsOperator not initialized");
    }

    std::vector<double> result;
    result.reserve(obs.size());

    // Get state dimensions
    size_t nx = state.geometry().x_dim();
    size_t ny = state.geometry().y_dim();

    // Apply observation operator to each observation point
    for (const auto& obs_point : obs) {
      if (!obs_point.is_valid) {
        result.push_back(0.0);  // Invalid observations get zero
        continue;
      }

      // Convert observation location to grid coordinates
      double x, y;
      if (obs_point.location.getCoordinateSystem() ==
          CoordinateSystem::GEOGRAPHIC) {
        auto [lat, lon, level] = obs_point.location.getGeographicCoords();
        x = lon;
        y = lat;
      } else {
        // For non-geographic coordinates, use grid coordinates directly
        auto [i, j] = obs_point.location.getGridCoords2D();
        x = static_cast<double>(i);
        y = static_cast<double>(j);
      }

      // Perform nearest-neighbor interpolation
      double interpolated_value = nearestNeighborInterpolation(state, x, y, nx, ny);
      result.push_back(interpolated_value);
    }

    return result;
  }

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

 private:
  /**
   * @brief Perform nearest-neighbor interpolation from grid to point
   *
   * @param state Model state
   * @param x Longitude coordinate
   * @param y Latitude coordinate
   * @param nx Number of grid points in x direction
   * @param ny Number of grid points in y direction
   * @return Value at the nearest grid point
   */
  double nearestNeighborInterpolation(const StateBackend& state, double x, double y,
                            size_t nx, size_t ny) const {
    // Clamp coordinates to grid bounds
    x = std::max(0.0, std::min(static_cast<double>(nx - 1), x));
    y = std::max(0.0, std::min(static_cast<double>(ny - 1), y));

    // Find nearest grid point (round to nearest integer)
    size_t i = static_cast<size_t>(std::round(x));
    size_t j = static_cast<size_t>(std::round(y));

    // Ensure indices are within bounds
    if (i >= nx) i = nx - 1;
    if (j >= ny) j = ny - 1;

    // Return exact value at the grid point using grid coordinates
    return state.at(std::make_pair(static_cast<int>(i), static_cast<int>(j)));
  }

  bool initialized_ = false;                      ///< Initialization status
  std::vector<std::string> required_state_vars_;  ///< Required state variables
  std::vector<std::string>
      required_obs_vars_;  ///< Required observation variables
};

}  // namespace metada::backends::common::obsoperator 