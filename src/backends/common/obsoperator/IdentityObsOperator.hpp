/**
 * @file IdentityObsOperator.hpp
 * @brief Identity observation operator backend implementation
 * @ingroup backends
 * @author Metada Framework Team
 *
 * @details
 * This class provides an identity observation operator backend interface that maps 
 * any state backend to any observation backend. The observation operator performs nearest-neighbor 
 * interpolation from the model grid to observation locations.
 */

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "Location.hpp"
#include "PointObservation.hpp"

namespace metada::backends::common::obsoperator {

using framework::CoordinateSystem;

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
 * - Works with any state and observation backends that provide required interfaces
 *
 * @tparam StateBackend Type of state backend to operate on
 * @tparam ObsBackend Type of observation backend to operate on
 */
template <typename StateBackend, typename ObsBackend>
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
                            const ObsBackend& obs) const {
    if (!isInitialized()) {
      throw std::runtime_error("IdentityObsOperator not initialized");
    }

    std::vector<double> result;
    result.reserve(obs.size());

    // Get state dimensions
    size_t nx = state.geometry().x_dim();
    size_t ny = state.geometry().y_dim();
    size_t nz = state.geometry().z_dim();

    // Apply observation operator to each observation point
    for (const auto& obs_point : obs) {
      if (!obs_point.is_valid) {
        result.push_back(0.0);  // Invalid observations get zero
        continue;
      }

      // Convert observation location to grid coordinates
      double x, y, z;
      if (obs_point.location.getCoordinateSystem() ==
          CoordinateSystem::GEOGRAPHIC) {
        auto [lat, lon, level] = obs_point.location.getGeographicCoords();
        
        // Convert geographic coordinates to grid coordinates
        auto [grid_i, grid_j, grid_k] = convertGeographicToGrid(lat, lon, level, state.geometry());
        x = static_cast<double>(grid_i);
        y = static_cast<double>(grid_j);
        z = static_cast<double>(grid_k);
      } else {
        // For non-geographic coordinates, use grid coordinates directly
        auto [i, j, k] = obs_point.location.getGridCoords();
        x = static_cast<double>(i);
        y = static_cast<double>(j);
        z = static_cast<double>(k);
      }

      // Perform nearest-neighbor interpolation
      double interpolated_value = nearestNeighborInterpolation(state, x, y, z, nx, ny, nz);
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
   * @brief Convert geographic coordinates to grid coordinates
   *
   * @details This method performs a nearest-neighbor search to find the closest
   * grid point to the given geographic coordinates. It requires the geometry
   * backend to provide coordinate lookup capabilities.
   *
   * @param lat Latitude in degrees
   * @param lon Longitude in degrees
   * @param level Vertical level/pressure
   * @param geometry Reference to the geometry backend
   * @return Tuple of (i, j, k) grid coordinates
   */
  template <typename GeometryBackend>
  std::tuple<size_t, size_t, size_t> convertGeographicToGrid(double lat, double lon, 
                                                   double level,
                                                   const GeometryBackend& geometry) const {
    // Try to use geometry-specific coordinate conversion if available
    if constexpr (requires { geometry.unstaggered_info(); }) {
      // For WRF-type geometry with coordinate arrays
      const auto& grid_info = geometry.unstaggered_info();
      
      if (grid_info.has_2d_coords()) {
        // Find closest grid point by searching through coordinate arrays
        double min_dist = std::numeric_limits<double>::max();
        size_t min_idx = 0;

        for (size_t idx = 0; idx < grid_info.longitude_2d.size(); ++idx) {
          double grid_lon = grid_info.longitude_2d[idx];
          double grid_lat = grid_info.latitude_2d[idx];

          // Calculate Euclidean distance in lat/lon space
          double dist = std::sqrt((lon - grid_lon) * (lon - grid_lon) +
                                 (lat - grid_lat) * (lat - grid_lat));

          if (dist < min_dist) {
            min_dist = dist;
            min_idx = idx;
          }
        }

        // Convert linear index to 2D grid coordinates
        size_t j = min_idx / grid_info.nx;
        size_t i = min_idx % grid_info.nx;
        
        // Find closest vertical level
        size_t k = 0;
        if (grid_info.has_vertical_coords()) {
          double min_vert_dist = std::numeric_limits<double>::max();
          for (size_t z = 0; z < grid_info.vertical_coords.size(); ++z) {
            double dist = std::abs(level - grid_info.vertical_coords[z]);
            if (dist < min_vert_dist) {
              min_vert_dist = dist;
              k = z;
            }
          }
        }
        
        return std::make_tuple(i, j, k);
      }
    }
    
    // Fallback: Search through geometry locations (slower but more general)
    double min_dist = std::numeric_limits<double>::max();
    size_t best_i = 0, best_j = 0, best_k = 0;
    
    for (size_t k = 0; k < geometry.z_dim(); ++k) {
      for (size_t j = 0; j < geometry.y_dim(); ++j) {
        for (size_t i = 0; i < geometry.x_dim(); ++i) {
          try {
            // Try to get geographic coordinates for this grid point
            auto grid_location = geometry.getLocation(i, j, k);
            if (grid_location.getCoordinateSystem() == CoordinateSystem::GEOGRAPHIC) {
              auto [grid_lat, grid_lon, grid_level] = grid_location.getGeographicCoords();
              
              // Calculate distance in 3D space (lat, lon, level)
              double horizontal_dist = std::sqrt((lon - grid_lon) * (lon - grid_lon) +
                                               (lat - grid_lat) * (lat - grid_lat));
              double vertical_dist = std::abs(level - grid_level);
              
              // Combined distance (you might want to weight these differently)
              double dist = horizontal_dist + vertical_dist * 0.001; // Scale vertical distance
              
              if (dist < min_dist) {
                min_dist = dist;
                best_i = i;
                best_j = j;
                best_k = k;
              }
            }
          } catch (...) {
            // Skip points where geographic coordinates are not available
            continue;
          }
        }
      }
    }
    
    return std::make_tuple(best_i, best_j, best_k);
  }

  /**
   * @brief Perform nearest-neighbor interpolation from grid to point
   *
   * @param state Model state
   * @param x Grid i-coordinate (already converted from geographic if needed)
   * @param y Grid j-coordinate (already converted from geographic if needed)
   * @param z Grid k-coordinate (already converted from geographic if needed)
   * @param nx Number of grid points in x direction
   * @param ny Number of grid points in y direction
   * @param nz Number of grid points in z direction
   * @return Value at the nearest grid point
   */
  double nearestNeighborInterpolation(const StateBackend& state, double x, double y, double z,
                            size_t nx, size_t ny, size_t nz) const {
    // Clamp coordinates to grid bounds
    x = std::max(0.0, std::min(static_cast<double>(nx - 1), x));
    y = std::max(0.0, std::min(static_cast<double>(ny - 1), y));
    z = std::max(0.0, std::min(static_cast<double>(nz - 1), z));

    // Find nearest grid point (round to nearest integer)
    size_t i = static_cast<size_t>(std::round(x));
    size_t j = static_cast<size_t>(std::round(y));
    size_t k = static_cast<size_t>(std::round(z));

    // Ensure indices are within bounds
    if (i >= nx) i = nx - 1;
    if (j >= ny) j = ny - 1;
    if (k >= nz) k = nz - 1;

    // Convert 3D grid coordinates to linear index using row-major order
    // linear_index = k * (nx * ny) + j * nx + i
    size_t linear_index = k * (nx * ny) + j * nx + i;

    // Use linear indexing to avoid coordinate system mismatch
    // This works with all backends that support operator[] or at(size_t)
    return state[linear_index];
  }

  bool initialized_ = false;                      ///< Initialization status
  std::vector<std::string> required_state_vars_;  ///< Required state variables
  std::vector<std::string>
      required_obs_vars_;  ///< Required observation variables
};

}  // namespace metada::backends::common::obsoperator 